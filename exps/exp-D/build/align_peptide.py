#!/usr/bin/env python3
"""Kabsch-align tleap-built peptide to 3IOL coordinates."""
import numpy as np
import parmed as pmd

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
THREEIOL = "/home/scroll/personal/semaglutide-in-silico/exps/exp-B/structures/3IOL.pdb"

# Load reference peptide from 3IOL (chain B only)
ref = pmd.load_file(THREEIOL)
ref_pep = ref[[a for a in ref.atoms if a.residue.chain == 'B']]
print(f"Reference peptide: {len(ref_pep.atoms)} atoms, "
      f"{len(set(r.idx for r in ref_pep.residues))} residues")

# Build tleap peptide sequence
# GLP-1(10-35) from 3IOL: GTFTSDVSSYLEGQAAKEFIAWLVKGRG... wait, need to check
# Chain B residues from 3IOL:
# GLY10,THR11,PHE12,THR13,SER14,ASP15,VAL16,SER17,SER18,TYR19,
# LEU20,GLU21,GLY22,GLN23,ALA24,ALA25,LYS26,GLU27,PHE28,ILE29,
# ALA30,TRP31,LEU32,VAL33,LYS34,GLY35
# After Lys34→Arg mutation: ...GLY35 with ARG at position 34
seq_3l = ['GLY','THR','PHE','THR','SER','ASP','VAL','SER','SER','TYR',
          'LEU','GLU','GLY','GLN','ALA','ALA','LYS','GLU','PHE','ILE',
          'ALA','TRP','LEU','VAL','ARG','GLY']

# Remove waters from reference
ref_pep = ref_pep[[a for a in ref_pep.atoms if a.residue.name not in ('HOH','WAT')]]
print(f"Ref peptide (no water): {len(ref_pep.atoms)} atoms")

# Get reference backbone coords (N, CA, C)
ref_bb = []
for r in ref_pep.residues:
    for a in r.atoms:
        if a.name in ('N', 'CA', 'C'):
            ref_bb.append(a)
print(f"Ref backbone atoms: {len(ref_bb)}")

# Generate tleap input that builds the peptide from sequence
leap_in = f"{TLEAP}/build_pep.in"
with open(leap_in, 'w') as f:
    f.write("""source leaprc.protein.ff14SB
pep = sequence { GLY THR PHE THR SER ASP VAL SER SER TYR LEU GLU GLY GLN ALA ALA LYS GLU PHE ILE ALA TRP LEU VAL ARG GLY }
savePdb pep pep_tleap.pdb
quit
""")

print(f"Wrote: {leap_in}")
print("Run: tleap -f build_pep.in")
print("Then run align script step 2 to Kabsch align pep_tleap.pdb to 3IOL")
