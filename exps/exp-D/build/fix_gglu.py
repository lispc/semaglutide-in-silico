#!/usr/bin/env python3
"""Rebuild gglu with rotated C18 tail to avoid ECD clash."""
import numpy as np, os, subprocess as sp

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"
os.chdir(TLEAP)

NZ_POS = np.array([-7.571, -0.448, 1.906])
TARGET_DIR = np.array([0.753, -0.615, 0.233])
DELETE = {'C1','C2','O3','N4','C5','C6','C7','C8','C9','C32','O33','N34','C35'}
RENAME = {'N10':'N', 'C11':'C', 'O12':'O'}

# Parse bcc mol2
atoms = {}
in_atom = False
with open(f"{BUILD}/gglu_bcc.mol2") as f:
    for line in f:
        ls = line.strip()
        if '@<TRIPOS>ATOM' in ls: in_atom = True; continue
        if in_atom and ls.startswith('@'): break
        if in_atom and ls:
            p = ls.split()
            if len(p) >= 9:
                atoms[int(p[0])] = {'name': p[1], 'x': float(p[2]), 'y': float(p[3]),
                                    'z': float(p[4]), 'charge': float(p[8])}

# Identify atoms to keep
keep_ids = sorted([aid for aid in atoms if atoms[aid]['name'] not in DELETE])

# Assign types
for aid in keep_ids:
    old_name = atoms[aid]['name']
    if old_name in RENAME:
        atoms[aid]['name'] = RENAME[old_name]
    name = atoms[aid]['name']; elem = name[0]
    if name == 'C': atoms[aid]['atype'] = 'c'
    elif name == 'O': atoms[aid]['atype'] = 'o'
    elif elem == 'N': atoms[aid]['atype'] = 'n'
    elif elem == 'H': atoms[aid]['atype'] = 'hc'
    else: atoms[aid]['atype'] = 'c3'

# Position: align N->C with TARGET_DIR
n_id, c_id = keep_ids[0], keep_ids[1]
n_pos = np.array([atoms[n_id]['x'], atoms[n_id]['y'], atoms[n_id]['z']])
c_pos = np.array([atoms[c_id]['x'], atoms[c_id]['y'], atoms[c_id]['z']])
n_to_c = c_pos - n_pos

v1 = n_to_c / np.linalg.norm(n_to_c); v2 = TARGET_DIR
cos_a = np.dot(v1, v2)
if abs(cos_a-1) < 1e-8:
    R = np.eye(3)
else:
    k = np.cross(v1, v2); k /= np.linalg.norm(k)
    K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    R = np.eye(3) + np.sin(np.arccos(cos_a))*K + (1-cos_a)*K@K

for aid in keep_ids:
    pos = np.array([atoms[aid]['x'], atoms[aid]['y'], atoms[aid]['z']])
    new_pos = R @ (pos - n_pos) + NZ_POS
    atoms[aid]['x'], atoms[aid]['y'], atoms[aid]['z'] = new_pos

# Write mol2
natoms = len(keep_ids)
out_mol2 = "lnk_gglu_pos.mol2"
with open(out_mol2, 'w') as f:
    f.write("@<TRIPOS>MOLECULE\nLNK\n")
    f.write(f" {natoms:5d}     0     0     0     0\nSMALL\nGAFF2\nlinker-C18\n\n@<TRIPOS>ATOM\n")
    for new_id, old_id in enumerate(keep_ids, 1):
        a = atoms[old_id]
        f.write(f"{new_id:6d} {a['name']:6s} {a['x']:10.4f} {a['y']:10.4f} "
                f"{a['z']:10.4f} {a['atype']:4s} 1 LNK {a['charge']:.6f}\n")
    f.write("@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n     1 LNK         1 TEMP              0 ****  ****    0 ROOT\n")

new_n = np.array([atoms[n_id]['x'], atoms[n_id]['y'], atoms[n_id]['z']])
new_c = np.array([atoms[c_id]['x'], atoms[c_id]['y'], atoms[c_id]['z']])
d = np.linalg.norm(new_c - NZ_POS)
print(f"NZ->C={d:.2f}A, dot(target)={np.dot((new_c-NZ_POS)/d, TARGET_DIR):.3f}")
print(f"Saved {out_mol2} ({natoms} atoms)")

# Build tleap
tleap_input = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p
loadAmberParams lya_link.frcmod
complex = loadPdb ecd_pep.pdb
LNK = loadMol2 lnk_gglu_pos.mol2
sys = combine { complex LNK }
remove sys sys.117.HZ1
remove sys sys.117.HZ2
remove sys sys.117.HZ3
remove sys sys.128.N
remove sys sys.127
bond sys.117.NZ sys.128.C
solvateOct sys TIP3PBOX 10.0
addIonsRand sys Na+ 0
saveAmberParm sys gglu.prmtop gglu.inpcrd
quit
"""
result = sp.run(['tleap', '-f', '-'], input=tleap_input, capture_output=True,
                text=True, timeout=120, env={**os.environ,
                'PATH': '/home/scroll/miniforge3/envs/cgas-md/bin:' + os.environ.get('PATH','')})
for line in (result.stdout + result.stderr).split('\n'):
    if 'Error' in line or 'FATAL' in line or 'Exiting' in line:
        print(line)

size = os.path.getsize('gglu.prmtop') if os.path.exists('gglu.prmtop') else 0
print(f"gglu.prmtop: {size/1e6:.1f}MB")
