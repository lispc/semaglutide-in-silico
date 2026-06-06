#!/usr/bin/env python3
"""
Build minimal full semaglutide-GLP-1R complex for exp-F Phase 0 validation.

Strategy:
1. Extract receptor + Gs from 7KI0 (chains R, A, B, G, N; remove P chain and waters)
2. Extract full semaglutide (peptide + LNK) from exp-D fixed structure
3. Kabsch-align exp-D peptide backbone (residues 101-126, GLY10-GLY35)
   onto 7KI0 P-chain backbone (GLY10-GLY35)
4. Combine into single PDB
5. Solvate with OpenMM

This gives us: full GLP-1R (ECD+TMD) + Gs + nanobody + full semaglutide
in the correct cryo-EM pose, with the complete linker-C18 diacid from exp-D.
"""
import numpy as np
import os, sys

REPO = "/home/scroll/personal/semaglutide-in-silico"
KI0 = f"{REPO}/exps/exp-B/structures/7KI0.pdb"
EXPD = f"{REPO}/exps/exp-D/tleap/gglu_2oeg_fixed.pdb"
OUT = f"{REPO}/exps/exp-F/structures"
os.makedirs(OUT, exist_ok=True)


def parse_pdb_atoms(path):
    """Parse ATOM/HETATM records. Returns list of dicts."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            rec = line[0:6].strip()
            serial = int(line[6:11])
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append({
                "rec": rec, "serial": serial, "name": name,
                "resname": resname, "chain": chain, "resid": resid,
                "x": x, "y": y, "z": z, "line": line.rstrip("\n")
            })
    return atoms


def kabsch_align(mobile, target):
    """Kabsch alignment: returns R, t such that mobile @ R + t ≈ target.
    mobile, target: (N, 3) arrays."""
    assert mobile.shape == target.shape
    # Center
    c_mobile = mobile.mean(axis=0)
    c_target = target.mean(axis=0)
    P = mobile - c_mobile
    Q = target - c_target
    # Covariance
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    # Correct reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = c_target - c_mobile @ R
    return R, t


# === 1. Parse structures ===
print("Parsing 7KI0...")
atoms_7ki0 = parse_pdb_atoms(KI0)
print(f"  Total atoms: {len(atoms_7ki0)}")

print("Parsing exp-D fixed...")
atoms_expd = parse_pdb_atoms(EXPD)
print(f"  Total atoms: {len(atoms_expd)}")

# === 2. Extract receptor from 7KI0 (R, A, B, G, N chains; exclude P and waters) ===
receptor_atoms = [
    a for a in atoms_7ki0
    if a["chain"] in ("R", "A", "B", "G", "N")
    and a["resname"] not in ("HOH", "WAT")
]
print(f"Receptor atoms: {len(receptor_atoms)}")

# === 3. Extract full semaglutide from exp-D (residues 101-127, exclude waters) ===
ligand_atoms = [
    a for a in atoms_expd
    if 101 <= a["resid"] <= 127
    and a["resname"] not in ("HOH", "WAT", "Na+")
]
print(f"Ligand atoms from exp-D: {len(ligand_atoms)}")

# === 4. Extract backbone CA for alignment ===
# 7KI0 P chain: GLY 10 - GLY 35
ki0_ca = []
for a in atoms_7ki0:
    if a["chain"] == "P" and a["name"] == "CA" and 10 <= a["resid"] <= 35:
        ki0_ca.append((a["resid"], np.array([a["x"], a["y"], a["z"]])))
ki0_ca.sort(key=lambda x: x[0])
ki0_coords = np.array([c for _, c in ki0_ca])
print(f"7KI0 P-chain CA atoms for alignment: {len(ki0_ca)} (GLY10-GLY35)")

# exp-D: residues 101-126 correspond to GLY10-GLY35
expd_ca = []
for a in ligand_atoms:
    if a["name"] == "CA" and 101 <= a["resid"] <= 126:
        expd_ca.append((a["resid"], np.array([a["x"], a["y"], a["z"]])))
expd_ca.sort(key=lambda x: x[0])
expd_coords = np.array([c for _, c in expd_ca])
print(f"exp-D ligand CA atoms for alignment: {len(expd_ca)} (res 101-126)")

assert len(ki0_ca) == len(expd_ca), f"Mismatch: {len(ki0_ca)} vs {len(expd_ca)}"

# === 5. Kabsch alignment ===
R, t = kabsch_align(expd_coords, ki0_coords)
rmsd_before = np.sqrt(((expd_coords - ki0_coords) ** 2).sum() / len(expd_coords))
print(f"RMSD before alignment: {rmsd_before:.2f} Å")

# Apply to all ligand atoms
for a in ligand_atoms:
    v = np.array([a["x"], a["y"], a["z"]])
    v_new = v @ R + t
    a["x"], a["y"], a["z"] = v_new

# Verify aligned CA RMSD
aligned_ca = []
for a in ligand_atoms:
    if a["name"] == "CA" and 101 <= a["resid"] <= 126:
        aligned_ca.append(np.array([a["x"], a["y"], a["z"]]))
aligned_ca = np.array(aligned_ca)
rmsd_after = np.sqrt(((aligned_ca - ki0_coords) ** 2).sum() / len(ki0_coords))
print(f"RMSD after alignment: {rmsd_after:.2f} Å")

# === 6. Write combined PDB ===
# Re-serialize receptor atoms, then ligand atoms with new coords
output_path = f"{OUT}/7ki0_receptor_plus_expD_ligand.pdb"
with open(output_path, "w") as f:
    f.write("REMARK  Combined 7KI0 receptor + exp-D full semaglutide\n")
    f.write("REMARK  Receptor: 7KI0 chains R,A,B,G,N (GLP-1R + Gs + nanobody)\n")
    f.write("REMARK  Ligand: exp-D gglu_2oeg_fixed residues 101-127 (GLP-1 10-35 + LNK)\n")
    f.write("REMARK  Aligned by Kabsch on GLY10-GLY35 backbone CA\n")
    f.write(f"REMARK  Alignment RMSD: {rmsd_after:.2f} Å\n")
    
    serial = 1
    # Receptor
    for a in receptor_atoms:
        f.write(f"ATOM  {serial:5d} {a['name']:4s} {a['resname']:3s} {a['chain']:1s}{a['resid']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['name'][0]:1s}  \n")
        serial += 1
    
    # Ligand (chain P)
    for a in ligand_atoms:
        f.write(f"ATOM  {serial:5d} {a['name']:4s} {a['resname']:3s} P{a['resid']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['name'][0]:1s}  \n")
        serial += 1
    
    f.write("END\n")

print(f"\nWritten: {output_path}")
print(f"  Receptor atoms: {len(receptor_atoms)}")
print(f"  Ligand atoms: {len(ligand_atoms)}")
print(f"  Total: {len(receptor_atoms) + len(ligand_atoms)}")

# === 7. Quick check: NZ-C distance in ligand ===
nz = None
c0 = None
for a in ligand_atoms:
    if a["resname"] == "LYS" and a["name"] == "NZ":
        nz = np.array([a["x"], a["y"], a["z"]])
    if a["resname"] == "LNK" and a["name"] == "C":
        c0 = np.array([a["x"], a["y"], a["z"]])
if nz is not None and c0 is not None:
    dist = np.linalg.norm(nz - c0)
    print(f"  NZ-C (Lys26-LNK amide) distance: {dist:.2f} Å (expected ~1.5 Å)")
