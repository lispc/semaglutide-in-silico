#!/usr/bin/env python3
"""
Build ECD minimal model v2:
- ECD truncated to residue 128 (true ECD boundary)
- Peptide from 4ZGM crystal structure (more accurate than 7OR0)
- LNK translated to connect to 4ZGM Lys26
- HSA placed at +120 Å in Z
"""
import numpy as np

def read_pdb_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'record': line[:6].strip(),
                    'atomnum': int(line[6:11]),
                    'atomname': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resnum': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[77:78].strip() if len(line) > 77 else '',
                }
                atoms.append(atom)
    return atoms

def kabsch_align(mobile, target):
    """Kabsch alignment. Returns (transformed_coords, R, mobile_cent, target_cent)"""
    mobile_cent = np.mean(mobile, axis=0)
    target_cent = np.mean(target, axis=0)
    mobile_c = mobile - mobile_cent
    target_c = target - target_cent
    H = mobile_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    transformed = (mobile_c @ R) + target_cent
    return transformed, R, mobile_cent, target_cent

# === 1. Load 4ZGM structures ===
print("[1/6] Loading 4ZGM ECD + peptide...")
zgm_ecd = read_pdb_atoms('/home/scroll/personal/semaglutide-in-silico/exps/exp-B/structures/ecd_clean.pdb')
zgm_pep = read_pdb_atoms('/home/scroll/personal/semaglutide-in-silico/exps/exp-B/structures/pep_4zgm_clean.pdb')
print(f"  4ZGM ECD: {len(zgm_ecd)} atoms, res {min(a['resnum'] for a in zgm_ecd)}-{max(a['resnum'] for a in zgm_ecd)}")
print(f"  4ZGM peptide: {len(zgm_pep)} atoms, res {min(a['resnum'] for a in zgm_pep)}-{max(a['resnum'] for a in zgm_pep)}")

# === 2. Load 7OR0 ECD ===
print("[2/6] Loading 7OR0 ECD...")
or0_all = read_pdb_atoms('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/membrane_build/protein_final.pdb')
or0_ecd = [a for a in or0_all if a['chain'] == 'R' and 29 <= a['resnum'] <= 128]
print(f"  7OR0 ECD: {len(or0_ecd)} atoms, res {min(a['resnum'] for a in or0_ecd)}-{max(a['resnum'] for a in or0_ecd)}")

# === 3. Align 4ZGM ECD to 7OR0 ECD using CA atoms ===
print("[3/6] Aligning 4ZGM to 7OR0...")
zgm_ca = [(a['resnum'], np.array([a['x'], a['y'], a['z']])) for a in zgm_ecd if a['atomname'] == 'CA']
or0_ca = [(a['resnum'], np.array([a['x'], a['y'], a['z']])) for a in or0_ecd if a['atomname'] == 'CA']

common_res = sorted(set(r for r, _ in zgm_ca) & set(r for r, _ in or0_ca))
print(f"  Common CA atoms: {len(common_res)}")

zgm_coords = np.array([c for r, c in zgm_ca if r in common_res])
or0_coords = np.array([c for r, c in or0_ca if r in common_res])

transformed, R, zgm_cent, or0_cent = kabsch_align(zgm_coords, or0_coords)
rmsd = np.sqrt(np.mean(np.sum((transformed - or0_coords)**2, axis=1)))
print(f"  CA RMSD after alignment: {rmsd:.3f} Å")

# Apply transformation to all 4ZGM atoms (ECD + peptide)
def transform_atom(a):
    p = np.array([a['x'], a['y'], a['z']])
    p_new = (p - zgm_cent) @ R + or0_cent
    a['x'], a['y'], a['z'] = p_new
    return a

for a in zgm_ecd:
    transform_atom(a)
for a in zgm_pep:
    transform_atom(a)

# Change chain IDs
for a in zgm_ecd:
    a['chain'] = 'R'
for a in zgm_pep:
    a['chain'] = 'P'

# === 4. Find Lys26 NZ in transformed 4ZGM peptide ===
lys26_nz = None
for a in zgm_pep:
    if a['resnum'] == 26 and a['atomname'] == 'NZ':
        lys26_nz = np.array([a['x'], a['y'], a['z']])
        break
if lys26_nz is None:
    raise ValueError("Lys26 NZ not found in 4ZGM peptide!")
print(f"  Transformed Lys26 NZ: {lys26_nz}")

# === 5. Load and position LNK from 7OR0 ===
print("[4/6] Positioning LNK...")
or0_lnk = [a for a in or0_all if a['chain'] == 'P' and a['resnum'] == 127]
print(f"  7OR0 LNK: {len(or0_lnk)} atoms")

# Find 7OR0 Lys117 NZ (Lys26 equivalent in 7OR0 numbering)
or0_lys117_nz = None
for a in or0_all:
    if a['chain'] == 'P' and a['resnum'] == 117 and a['atomname'] == 'NZ':
        or0_lys117_nz = np.array([a['x'], a['y'], a['z']])
        break
if or0_lys117_nz is None:
    raise ValueError("Lys117 NZ not found in 7OR0!")

# Find LNK C atom (first atom in mol2)
lnc = None
for a in or0_lnk:
    if a['atomname'] == 'C':
        lnc = np.array([a['x'], a['y'], a['z']])
        break
if lnc is None:
    raise ValueError("LNK C not found!")

# Vector from 7OR0 Lys117 NZ to LNK C
vector = lnc - or0_lys117_nz
print(f"  7OR0 Lys117 NZ -> LNK C vector: {vector}, length: {np.linalg.norm(vector):.3f} Å")

# Translate LNK so that its C is at the same relative position to 4ZGM Lys26 NZ
target_lnc = lys26_nz + vector
translation = target_lnc - lnc

for a in or0_lnk:
    p = np.array([a['x'], a['y'], a['z']])
    p_new = p + translation
    a['x'], a['y'], a['z'] = p_new
    a['chain'] = 'L'
    a['resnum'] = 129

# Verify distance
new_lnc = lnc + translation
print(f"  New LNK C position: {new_lnc}")
print(f"  Distance Lys26 NZ -> LNK C: {np.linalg.norm(lys26_nz - new_lnc):.3f} Å")

# === 6. Load and translate HSA ===
print("[5/6] Loading HSA...")
hsa_atoms = read_pdb_atoms('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/minimal_model/hsa_clean_no_myr.pdb')
print(f"  HSA: {len(hsa_atoms)} atoms")

offset = np.array([0.0, 0.0, 120.0])
for a in hsa_atoms:
    a['x'] += offset[0]
    a['y'] += offset[1]
    a['z'] += offset[2]

# === 7. Merge all ===
print("[6/6] Merging and writing...")
all_atoms = zgm_ecd + zgm_pep + or0_lnk + hsa_atoms

out_path = '/home/scroll/personal/semaglutide-in-silico/exps/exp-F/minimal_model/complex_ecd_v2.pdb'
with open(out_path, 'w') as f:
    serial = 1
    for a in all_atoms:
        # Skip hydrogen atoms - let tleap add them with correct names
        if a['element'] == 'H':
            continue
        
        # Ensure atom name is 4 chars, left-aligned for 3-4 char names, right-aligned for 1-2 char
        name = a['atomname']
        if len(name) <= 2:
            name_str = f" {name:>2s} "
        else:
            name_str = f"{name:<4s}"
        
        f.write(f"ATOM  {serial:5d} {name_str} {a['resname']:3s} {a['chain']:1s}{a['resnum']:4d}    {a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00           {a['element']:>1s}\n")
        serial += 1
    
    # TER records between chains
    # After ECD
    f.write(f"TER   {serial:5d}      {zgm_ecd[-1]['resname']:3s} {zgm_ecd[-1]['chain']:1s}{zgm_ecd[-1]['resnum']:4d}\n")
    serial += 1
    # After peptide (no TER before LNK since they should be connected)
    # After LNK
    f.write(f"TER   {serial:5d}      {or0_lnk[-1]['resname']:3s} {or0_lnk[-1]['chain']:1s}{or0_lnk[-1]['resnum']:4d}\n")
    serial += 1
    # After HSA
    f.write(f"TER   {serial:5d}      {hsa_atoms[-1]['resname']:3s} {hsa_atoms[-1]['chain']:1s}{hsa_atoms[-1]['resnum']:4d}\n")
    serial += 1

print(f"  Total atoms written: {serial-1}")
print(f"  Saved: {out_path}")

# Print summary
print("\n=== Summary ===")
print(f"ECD (chain R): {len(zgm_ecd)} atoms, res 29-128")
print(f"Peptide (chain P): {len(zgm_pep)} atoms, res 10-37 -> will renumber to 101-128")
print(f"LNK (chain L): {len(or0_lnk)} atoms, res 129")
print(f"HSA (chain A): {len(hsa_atoms)} atoms")
