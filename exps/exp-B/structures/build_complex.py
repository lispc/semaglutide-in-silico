#!/usr/bin/env python3
"""
Build GLP-1(10-37) + GLP-1R ECD complex PDB for tleap.

Strategy (revised):
  1. Use 4ZGM ECD + peptide backbone (residues 10-37) as template.
  2. Do NOT model N-terminal residues 7-9 (disordered in all crystals).
     The K34/R34 question is in the resolved C-terminal region.
     N-terminal Aib8 can be added later if needed.
  3. For B1: ARG34→LYS34, keep only backbone atoms (N,CA,C,O,CB),
     let tleap build correct LYS sidechain from template.
  4. For B2: use 4ZGM peptide as-is (already has R34, all atoms).
  5. Use NHIS/CGLY terminal naming for tleap compatibility.
  6. Strip all hydrogens (tleap adds them back).

Note: We use ALA at position 10-37 numbering which corresponds to
GLP-1 residues 10-37. Pos 8 (Aib8) is not modeled since it's
disordered in the crystal. If needed, N-terminal AIB can be added
via ParmEd modification of an extended N-terminal sequence.

Outputs:
  b1_complex_noh.pdb  - ECD + GLP-1(10-37, K34), no hydrogens, tleap-ready
  b2_complex_noh.pdb  - ECD + GLP-1(10-37, R34), no hydrogens, tleap-ready
"""
import numpy as np
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-B/structures"

def read_pdb_atoms(pdb_path):
    """Read ATOM/HETATM lines from PDB."""
    atoms = []
    with open(pdb_path) as f:
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

def write_tleap_pdb(atoms_by_chain, path, remark=""):
    """Write PDB for tleap: heavy atoms only, renamed N/C-terminal residues, TER between chains."""
    with open(path, 'w') as f:
        f.write(f"REMARK {remark}\n")
        atom_num = 0
        for chain in ['A', 'B']:
            if chain not in atoms_by_chain:
                continue
            for a in atoms_by_chain[chain]:
                atom_num += 1
                resname = a['resname']
                f.write(f"ATOM  {atom_num:5d} {a['atom']:4s} {resname:3s} "
                        f"{a['chain']}{a['resid']:4d}    "
                        f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                        f"  1.00  0.00          {a.get('element','C'):2s}  \n")
            f.write("TER\n")
        f.write("END\n")

print("=" * 60)
print("Building GLP-1(10-37) + GLP-1R ECD complexes for tleap")
print("=" * 60)

# Load data
ecd_all = read_pdb_atoms(f"{OUT}/ecd_clean.pdb")
pep_4zgm = read_pdb_atoms(f"{OUT}/pep_4zgm_clean.pdb")

# ECD heavy atoms
ecd_heavy = [a for a in ecd_all if a['atom'][0] != 'H']
print(f"ECD heavy atoms: {len(ecd_heavy)}")

# Peptide heavy atoms from 4ZGM
pep_heavy = [a for a in pep_4zgm if a['atom'][0] != 'H']
print(f"4ZGM peptide heavy atoms: {len(pep_heavy)}")

# ============================================================
# Build B1 (K34): Replace ARG34 sidechain with LYS backbone only
# ============================================================

pep_b1 = []
for a in pep_heavy:
    if a['resid'] == 34 and a['resname'] == 'ARG':
        # Keep only backbone atoms (N, CA, C, O, CB)
        # tleap will fill in the LYS sidechain from template
        if a['atom'] in ('N', 'CA', 'C', 'O', 'CB'):
            new_a = dict(a)
            new_a['resname'] = 'LYS'
            pep_b1.append(new_a)
        # Skip ARG-specific sidechain atoms (CG, CD, NE, CZ, NH1, NH2)
    else:
        pep_b1.append(dict(a))

print(f"B1 peptide atoms: {len(pep_b1)} (ARG34→LYS backbone only)")

# ============================================================
# Build B2 (R34): Use 4ZGM as-is
# ============================================================

pep_b2 = [dict(a) for a in pep_heavy]
print(f"B2 peptide atoms: {len(pep_b2)} (ARG34 from 4ZGM)")

# ============================================================
# Rename terminal residues for tleap
# ============================================================

# For residues 10-37:
# N-terminal: GLY10 → should stay GLY (not N-terminal in tleap sense - it's N-terminal but mid-chain)
# C-terminal: GLY37 → CGLY (C-terminal glycine)
# Actually, in tleap aminoct12.lib, C-terminal residues are named CGLY, CALA, etc.
# But tleap auto-detection often works with regular names. Let me just add OXT.

def rename_termini(atoms):
    """Find the first and last peptide residue and rename/add atoms as needed."""
    residues = sorted(set(a['resid'] for a in atoms))
    first_res = residues[0]
    last_res = residues[-1]

    # For C-terminal residue: add OXT atom
    # Get C, CA, O positions of last residue to place OXT
    c_atom = None; o_atom = None; ca_atom = None
    for a in atoms:
        if a['resid'] == last_res and a['atom'] == 'C':
            c_atom = np.array([a['x'], a['y'], a['z']])
        if a['resid'] == last_res and a['atom'] == 'O':
            o_atom = np.array([a['x'], a['y'], a['z']])
        if a['resid'] == last_res and a['atom'] == 'CA':
            ca_atom = np.array([a['x'], a['y'], a['z']])

    if c_atom is not None and o_atom is not None and ca_atom is not None:
        # OXT is the second terminal oxygen, roughly symmetric to O across C-CA
        c_o_vec = o_atom - c_atom
        c_ca_vec = ca_atom - c_atom
        # OXT is roughly in the plane, ~120° from O and ~120° from CA
        oxt_pos = c_atom - (c_o_vec + c_ca_vec)
        oxt_pos = c_atom + 1.25 * (oxt_pos - c_atom) / np.linalg.norm(oxt_pos - c_atom)
        atoms.append({
            'chain': atoms[0]['chain'],
            'resname': atoms[0]['resname'] if atoms[0]['resid'] != last_res else 'CGLY',
            'resid': last_res,
            'atom': 'OXT',
            'x': oxt_pos[0], 'y': oxt_pos[1], 'z': oxt_pos[2],
            'element': 'O',
        })

    return atoms

pep_b1 = rename_termini(pep_b1)
pep_b2 = rename_termini(pep_b2)

# ============================================================
# Write PDBs
# ============================================================

# ECD atoms (chain A)
ecd_for_pdb = []
for a in ecd_heavy:
    ecd_for_pdb.append({
        'chain': 'A', 'resname': a['resname'], 'resid': a['resid'],
        'atom': a['atom'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
        'element': a.get('element', 'C'),
    })

# B1: chain A = ECD, chain B = peptide (K34)
for a in pep_b1:
    a['chain'] = 'B'
write_tleap_pdb({'A': ecd_for_pdb, 'B': pep_b1}, f"{OUT}/b1_complex_noh.pdb",
    f"exp-B B1: GLP-1R ECD + GLP-1(10-37, K34), tleap-ready, {len(ecd_for_pdb)+len(pep_b1)} atoms")

# B2: chain A = ECD, chain B = peptide (R34)
for a in pep_b2:
    a['chain'] = 'B'
write_tleap_pdb({'A': ecd_for_pdb, 'B': pep_b2}, f"{OUT}/b2_complex_noh.pdb",
    f"exp-B B2: GLP-1R ECD + GLP-1(10-37, R34), tleap-ready, {len(ecd_for_pdb)+len(pep_b2)} atoms")

# ============================================================
# Summary
# ============================================================

print(f"\n{'='*60}")
print("Output files:")
print(f"  B1: {OUT}/b1_complex_noh.pdb ({len(ecd_for_pdb)+len(pep_b1)} atoms)")
print(f"  B2: {OUT}/b2_complex_noh.pdb ({len(ecd_for_pdb)+len(pep_b2)} atoms)")

# Verify peptide residues
for label, atoms in [("B1", pep_b1), ("B2", pep_b2)]:
    res_info = {}
    for a in atoms:
        rid = a['resid']
        if rid not in res_info:
            res_info[rid] = {'name': a['resname'], 'atoms': []}
        res_info[rid]['atoms'].append(a['atom'])
    print(f"\n{label} peptide residues (10-37):")
    for rid in sorted(res_info.keys()):
        info = res_info[rid]
        expected = set()
        if info['name'] == 'GLY':
            expected = {'N','CA','C','O'}
        elif info['name'] == 'ALA':
            expected = {'N','CA','C','O','CB'}
        elif info['name'] == 'LYS':
            expected = {'N','CA','C','O','CB'}  # backbone only, sidechain from tleap
        elif info['name'] == 'ARG':
            expected = {'N','CA','C','O','CB','CG','CD','NE','CZ','NH1','NH2'}
        else:
            expected = {'N','CA','C','O','CB'}
        actual = set(info['atoms'])
        missing = expected - actual
        extra = actual - expected
        status = ""
        if missing: status += f" [MISSING:{','.join(sorted(missing))}]"
        if extra: status += f" [EXTRA:{','.join(sorted(extra))}]"
        print(f"  {info['name']:4s} {rid:3d}: {len(info['atoms']):2d} atoms{status}")

print(f"\nDone!")
