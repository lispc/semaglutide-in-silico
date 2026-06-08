#!/usr/bin/env python3
import numpy as np

def modify_line(line, chain=None, x=None, y=None, z=None):
    line = line.rstrip('\n').ljust(80)
    chars = list(line)
    if chain is not None:
        chars[21] = chain
    if x is not None:
        chars[30:38] = list(f"{x:8.3f}")
    if y is not None:
        chars[38:46] = list(f"{y:8.3f}")
    if z is not None:
        chars[46:54] = list(f"{z:8.3f}")
    return ''.join(chars)

ecd_atoms = []
pep_atoms = []
lnk_atoms = []

with open('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/membrane_build/protein_final.pdb') as f:
    for line in f:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        chain = line[21]
        try:
            resnum = int(line[22:26])
        except:
            continue
        if chain == 'R' and 29 <= resnum <= 145:
            ecd_atoms.append(modify_line(line))
        elif chain == 'P' and resnum <= 126:
            pep_atoms.append(modify_line(line))
        elif chain == 'P' and resnum == 127:
            lnk_atoms.append(modify_line(line, chain='L'))

if not ecd_atoms:
    raise ValueError("No ECD atoms found!")
if not pep_atoms:
    raise ValueError("No peptide atoms found!")
if not lnk_atoms:
    raise ValueError("No LNK atoms found!")

hsa_atoms = []
offset = np.array([0.0, 0.0, 120.0])
with open('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/minimal_model/hsa_clean_no_myr.pdb') as f:
    for line in f:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        x = float(line[30:38]) + offset[0]
        y = float(line[38:46]) + offset[1]
        z = float(line[46:54]) + offset[2]
        hsa_atoms.append(modify_line(line, x=x, y=y, z=z))

if not hsa_atoms:
    raise ValueError("No HSA atoms found!")

with open('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/minimal_model/complex_ecd.pdb', 'w') as f:
    serial = [1]
    def write_chain(atoms, ter=True):
        for line in atoms:
            new_line = list(line.ljust(80))
            new_line[6:11] = list(f"{serial[0]:5d}")
            f.write(''.join(new_line) + '\n')
            serial[0] += 1
        if ter:
            last = atoms[-1]
            resName = last[17:20].strip()
            chain = last[21]
            resnum = int(last[22:26])
            f.write(f"TER   {serial[0]:5d}      {resName:3s} {chain:1s}{resnum:4d}\n")
            serial[0] += 1

    write_chain(ecd_atoms)
    write_chain(pep_atoms, ter=False)
    write_chain(lnk_atoms, ter=False)
    # TER after LNK
    last = lnk_atoms[-1]
    resName = last[17:20].strip()
    chain = last[21]
    resnum = int(last[22:26])
    f.write(f"TER   {serial[0]:5d}      {resName:3s} {chain:1s}{resnum:4d}\n")
    serial[0] += 1
    write_chain(hsa_atoms)

all_atoms = ecd_atoms + pep_atoms + lnk_atoms + hsa_atoms
zs = [float(line[46:54]) for line in all_atoms]
print(f"Total atoms: {len(all_atoms)}")
print(f"Z range: {min(zs):.1f} to {max(zs):.1f}, span={max(zs)-min(zs):.1f}")
