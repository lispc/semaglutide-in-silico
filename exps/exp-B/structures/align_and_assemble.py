#!/usr/bin/env python3
"""
Align tleap-built peptides to crystal structure coordinates and assemble
the full complex PDB for final tleap processing.

Steps:
  1. Read tleap-built peptide (extended chain from sequence command)
  2. Read crystal structure peptide coordinates (bound pose, residues 10-37)
  3. Kabsch-align tleap peptide residues 10-37 to crystal residues 10-37
  4. Apply transformation to all peptide atoms (including residues 7-9)
  5. Combine with ECD PDB
  6. Write final complex PDB ready for tleap solvation
"""
import numpy as np
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-B/structures"
TLEAP_DIR = f"{REPO}/exps/exp-B/tleap"

def read_pdb(path):
    """Read PDB ATOM/HETATM records. Returns list of dicts."""
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append({
                    'chain': line[21],
                    'resname': line[17:20].strip(),
                    'resid': int(line[22:26]),
                    'atom': line[12:16].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[76:78].strip() or 'C',
                })
    return atoms

def write_pdb(atoms_by_chain, path, remark=""):
    """Write PDB with TER between chains."""
    atom_num = 0
    with open(path, 'w') as f:
        f.write(f"REMARK {remark}\n")
        for chain in ['A', 'B']:
            if chain not in atoms_by_chain:
                continue
            for a in atoms_by_chain[chain]:
                atom_num += 1
                f.write(f"ATOM  {atom_num:5d} {a['atom']:4s} {a['resname']:3s} "
                        f"{a['chain']}{a['resid']:4d}    "
                        f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                        f"  1.00  0.00          {a.get('element','C'):2s}  \n")
            f.write("TER\n")
        f.write("END\n")

def kabsch(P, Q):
    """Find optimal rotation matrix to align P onto Q. P, Q are Nx3 arrays."""
    P_cent = P.mean(axis=0)
    Q_cent = Q.mean(axis=0)
    P_c = P - P_cent
    Q_c = Q - Q_cent
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Q_cent - R @ P_cent
    return R, t

print("=" * 60)
print("Aligning tleap peptides to crystal structure coordinates")
print("=" * 60)

# Load crystal structure coordinates for peptide residues 10-37
# These are extracted from 4ZGM pep_4zgm_clean.pdb
pep_crystal = read_pdb(f"{OUT}/pep_4zgm_clean.pdb")
# Heavy atoms only
pep_crystal_hv = [a for a in pep_crystal if a['atom'][0] != 'H']

# Map crystal residues to atom positions for Kabsch
# Crystal: residues 10-37, we use backbone N, CA, C for alignment
crystal_bb = {}  # (resid, atom) -> coords
for a in pep_crystal_hv:
    if a['atom'] in ('N', 'CA', 'C'):
        crystal_bb[(a['resid'], a['atom'])] = np.array([a['x'], a['y'], a['z']])

print(f"Crystal peptide backbone atoms: {len(crystal_bb)}")

# ============================================================
# Process B1 (K34) and B2 (R34)
# ============================================================

for label, tleap_pdb in [("b1", "pep_b1_tleap.pdb"), ("b2", "pep_b2_tleap.pdb")]:
    print(f"\n{'='*60}")
    print(f"Processing {label.upper()}")
    print(f"{'='*60}")

    # Read tleap-built peptide
    pep_tleap = read_pdb(f"{TLEAP_DIR}/{tleap_pdb}")
    print(f"  Tleap PDB: {len(pep_tleap)} atoms")

    # Tleap residues: 1=His7, 2=Ala8, ..., 4=Gly10, ..., 31=Gly37
    # Crystal residues: Gly10 at residue 10, etc.
    # Mapping: tleap_res = crystal_res - 6 (for residues 10-37 → 4-31)
    # Actually, in the tleap output PDB, the residues are numbered 1-31 sequentially.
    # Let me verify:
    tleap_resids = sorted(set(a['resid'] for a in pep_tleap))
    print(f"  Tleap residues: {tleap_resids}")
    tleap_sequences = {}
    for rid in tleap_resids:
        names = set(a['resname'] for a in pep_tleap if a['resid'] == rid)
        tleap_sequences[rid] = list(names)[0]
    print(f"  Sequence (first 5): {[(r, tleap_sequences[r]) for r in tleap_resids[:5]]}")
    print(f"  Sequence (last 5): {[(r, tleap_sequences[r]) for r in tleap_resids[-5:]]}")

    # Build Kabsch matrices
    # Tleap residue 4 = Gly10, tleap residue 28 = K34/R34, tleap residue 31 = Gly37
    # Crystal residue 10 = Gly10, crystal residue 34 = K34/R34, crystal residue 37 = Gly37
    offset = 10 - 4  # crystal_res = tleap_res + offset

    P = []  # tleap coords
    Q = []  # crystal coords
    for tleap_rid in range(4, 32):  # residues 4-31 = crystal 10-37
        crystal_rid = tleap_rid + 6  # 4→10, 31→37
        for atom in ('N', 'CA', 'C'):
            # Get tleap coord
            tleap_hits = [a for a in pep_tleap if a['resid'] == tleap_rid and a['atom'] == atom]
            crys_key = (crystal_rid, atom)
            if tleap_hits and crys_key in crystal_bb:
                P.append([tleap_hits[0]['x'], tleap_hits[0]['y'], tleap_hits[0]['z']])
                Q.append(crystal_bb[crys_key])

    P = np.array(P)
    Q = np.array(Q)
    print(f"  Alignment points: {len(P)}")

    if len(P) < 3:
        print("  ERROR: Not enough alignment points")
        continue

    # Compute RMSD before alignment
    rmsd_before = np.sqrt(((P - Q)**2).mean())
    print(f"  RMSD before alignment: {rmsd_before:.1f} Å")

    # Kabsch alignment
    R, t = kabsch(P, Q)
    rmsd_after = np.sqrt((( (P @ R.T + t) - Q)**2).mean())
    print(f"  RMSD after alignment: {rmsd_after:.3f} Å")

    # Apply transformation to all tleap peptide atoms
    pep_aligned = []
    for a in pep_tleap:
        if a['atom'][0] == 'H':
            continue  # skip hydrogens (tleap will add them)
        new_coord = np.array([a['x'], a['y'], a['z']]) @ R.T + t
        pep_aligned.append({
            'chain': 'B',
            'resname': a['resname'],
            'resid': a['resid'] + 6,  # renumber to GLP-1 7-37
            'atom': a['atom'],
            'x': new_coord[0], 'y': new_coord[1], 'z': new_coord[2],
            'element': a.get('element', 'C'),
        })

    # Verify: check distance between aligned peptide Gly10 CA and crystal Gly10 CA
    aligned_g10_ca = [a for a in pep_aligned if a['resid'] == 10 and a['atom'] == 'CA']
    if aligned_g10_ca and (10, 'CA') in crystal_bb:
        ca_aligned = np.array([aligned_g10_ca[0]['x'], aligned_g10_ca[0]['y'], aligned_g10_ca[0]['z']])
        ca_crystal = crystal_bb[(10, 'CA')]
        d = np.linalg.norm(ca_aligned - ca_crystal)
        print(f"  Gly10 CA alignment check: {d:.4f} Å (should be ~0)")

    # Load ECD heavy atoms
    ecd_all = read_pdb(f"{OUT}/ecd_clean.pdb")
    ecd_heavy = []
    for a in ecd_all:
        if a['atom'][0] != 'H':
            ecd_heavy.append({
                'chain': 'A',
                'resname': a['resname'],
                'resid': a['resid'],
                'atom': a['atom'],
                'x': a['x'], 'y': a['y'], 'z': a['z'],
                'element': a.get('element', 'C'),
            })

    # Check for clashes between peptide and ECD
    pep_coords = np.array([[a['x'], a['y'], a['z']] for a in pep_aligned])
    ecd_coords = np.array([[a['x'], a['y'], a['z']] for a in ecd_heavy])
    # Quick check: minimum distance between any peptide atom and any ECD atom
    # This is slow (N*M) but our system is small
    from scipy.spatial import cKDTree
    ecd_tree = cKDTree(ecd_coords)
    dists, _ = ecd_tree.query(pep_coords, k=1)
    close = (dists < 2.0).sum()
    clash = (dists < 1.0).sum()
    print(f"  Peptide-ECD contacts: {close} pairs < 2 Å, {clash} pairs < 1 Å")

    # Write combined complex PDB
    complex_atoms = {'A': ecd_heavy, 'B': pep_aligned}
    total = len(ecd_heavy) + len(pep_aligned)
    out_path = f"{OUT}/{label}_complex_noh.pdb"
    write_pdb(complex_atoms, out_path,
              f"exp-B {label.upper()}: GLP-1R ECD + GLP-1(7-37, ALA8) aligned, no H, {total} atoms")

    # Verify residues
    pep_resids = sorted(set(a['resid'] for a in pep_aligned))
    pep_first = pep_resids[0]
    pep_last = pep_resids[-1]
    pep_seq = []
    for rid in pep_resids:
        rnames = set(a['resname'] for a in pep_aligned if a['resid'] == rid)
        pep_seq.append((rid, list(rnames)[0]))
    print(f"  Peptide residues: {pep_first}-{pep_last}")
    print(f"  First 3: {pep_seq[:3]}, Last 3: {pep_seq[-3:]}")
    print(f"  Output: {out_path}")

print(f"\n{'='*60}")
print("Done! Ready for tleap complex solvation.")
print(f"Output: {OUT}/b1_complex_noh.pdb, {OUT}/b2_complex_noh.pdb")
