#!/usr/bin/env python3
"""
Build minimal functional model: HSA + semaglutide + receptor ECD.
No membrane, no Gs. ~170k atoms target.
"""
import sys, os, numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-F/minimal_model"
os.makedirs(OUT, exist_ok=True)

def load_pdb_atoms(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21] if len(line) > 21 else ' '
                resnum = int(line[22:26].strip()) if len(line) > 26 else 0
                atoms.append({
                    'record': line[:6].strip(),
                    'atomnum': int(line[6:11]),
                    'atomname': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': chain,
                    'resnum': resnum,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[77:78].strip() if len(line) > 77 else '',
                })
    return atoms

print("[1/5] Loading receptor + peptide from exp-F...")
prot_atoms = load_pdb_atoms(f"{REPO}/exps/exp-F/membrane_build/protein_final.pdb")

# Receptor: chains R,A,B,G,N; Peptide: chain P (res 101-126); LNK: chain P (res 127)
rec_atoms = [a for a in prot_atoms if a['chain'] in ('R', 'A', 'B', 'G', 'N')]
pep_atoms = [a for a in prot_atoms if a['chain'] == 'P' and a['resnum'] <= 126]
lnk_atoms = [a for a in prot_atoms if a['chain'] == 'P' and a['resnum'] == 127]

print(f"  Receptor: {len(rec_atoms)} atoms")
print(f"  Peptide: {len(pep_atoms)} atoms")
print(f"  LNK: {len(lnk_atoms)} atoms")

print("[2/5] Loading HSA from exp-C...")
hsa_atoms = load_pdb_atoms(f"{REPO}/exps/exp-C/tleap/hsa_clean.pdb")
print(f"  HSA: {len(hsa_atoms)} atoms")

# === 3. Place HSA ===
# LNK tail end (C55/C58/O56/O57/O59)
lntail = [a for a in lnk_atoms if a['atomname'] in ['C55', 'C58', 'O56', 'O57', 'O59']]
if not lntail:
    # Fallback: use last few atoms of LNK
    lntail = lnk_atoms[-10:]
lnk_cent = np.mean([[a['x'], a['y'], a['z']] for a in lntail], axis=0)

hsa_cent = np.mean([[a['x'], a['y'], a['z']] for a in hsa_atoms], axis=0)

# Place HSA such that distance to LNK tail is ~8 nm
offset = np.array([0, 0, 80.0])  # 80 Å offset in z

print(f"  LNK tail centroid: {lnk_cent}")
print(f"  HSA centroid (original): {hsa_cent}")

for a in hsa_atoms:
    a['x'] += offset[0] - hsa_cent[0] + lnk_cent[0]
    a['y'] += offset[1] - hsa_cent[1] + lnk_cent[1]
    a['z'] += offset[2] - hsa_cent[2] + lnk_cent[2]

new_hsa_cent = np.mean([[a['x'], a['y'], a['z']] for a in hsa_atoms], axis=0)
print(f"  HSA centroid (translated): {new_hsa_cent}")
print(f"  Distance LNK-HSA: {np.linalg.norm(lnk_cent - new_hsa_cent):.1f} Å")

# === 4. Merge ===
print("[3/5] Merging structures...")

def write_pdb(atoms, chain_id, start_atomnum=1, start_resnum=1):
    lines = []
    atom_counter = start_atomnum
    res_map = {}
    for a in atoms:
        key = (a['resnum'], a['chain'])
        if key not in res_map:
            res_map[key] = start_resnum + len(res_map)
        new_resnum = res_map[key]
        
        atomnum_str = f"{atom_counter:>5d}" if atom_counter <= 99999 else f"{atom_counter:>5X}"
        resnum_str = f"{new_resnum:>4d}" if new_resnum <= 9999 else f"{new_resnum:>4X}"
        
        lines.append(f"ATOM  {atomnum_str} {a['atomname']:4s} {a['resname']:>3s} {chain_id:>1s}{resnum_str}    {a['x']:>8.3f}{a['y']:>8.3f}{a['z']:>8.3f}  1.00  0.00           {a['element']:>1s}\n")
        atom_counter += 1
    return lines, atom_counter

all_lines = []
atomnum = 1
resnum = 1

# Receptor (chain R)
lines, atomnum = write_pdb(rec_atoms, 'R', atomnum, resnum)
all_lines.extend(lines)
resnum += len(set((a['resnum'], a['chain']) for a in rec_atoms))

# Peptide (chain P)
lines, atomnum = write_pdb(pep_atoms, 'P', atomnum, resnum)
all_lines.extend(lines)
resnum += len(set((a['resnum'], a['chain']) for a in pep_atoms))

# LNK (chain L)
lines, atomnum = write_pdb(lnk_atoms, 'L', atomnum, resnum)
all_lines.extend(lines)
resnum += len(set((a['resnum'], a['chain']) for a in lnk_atoms))

# HSA (chain H)
lines, atomnum = write_pdb(hsa_atoms, 'H', atomnum, resnum)
all_lines.extend(lines)

all_lines.append("END\n")

with open(f"{OUT}/complex_minimal.pdb", 'w') as f:
    f.writelines(all_lines)

print(f"  Total atoms: {atomnum-1}")
print(f"  Saved: {OUT}/complex_minimal.pdb")

# === 5. tleap input ===
print("[4/5] Generating tleap input...")
lnk_mol2 = f"{REPO}/exps/exp-F/build/lnk_noh_zero.mol2"
frcmod = f"{REPO}/exps/exp-F/build/lya_link_final.frcmod"

tleap_in = f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadAmberParams {frcmod}
LNK = loadMol2 {lnk_mol2}

COMPLEX = loadPdb {OUT}/complex_minimal.pdb

SYS = combine {{COMPLEX}}

# Form amide bond between LYS P117 and LNK
# In the merged PDB, LYS is in chain P, LNK is in chain L
# Need to find correct tleap indices after loading
# We will inspect and fix manually if needed

desc SYS

solvateOct SYS TIP3PBOX 10.0
addIonsRand SYS Na+ 0
addIonsRand SYS Cl- 0

savePdb SYS {OUT}/system_minimal.pdb
saveAmberParm SYS {OUT}/system_minimal.prmtop {OUT}/system_minimal.inpcrd
quit
"""

with open(f"{OUT}/build_minimal.in", 'w') as f:
    f.write(tleap_in)

print(f"  Saved: {OUT}/build_minimal.in")
print("[5/5] Ready for tleap.")
print(f"  Run: cd {OUT} && tleap -f build_minimal.in")
