#!/usr/bin/env python3
"""Quick structural analysis of equilibrated exp-F system."""
import sys, math
import parmed as pmd
import numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
PRMTOP = f"{REPO}/exps/exp-F/tleap/system.prmtop"
EQUIL = f"{REPO}/exps/exp-F/md/rep1/equilibrated.pdb"
MIN = f"{REPO}/exps/exp-F/md/rep1/minimized.pdb"
REF = f"{REPO}/exps/exp-B/structures/7KI0.pdb"

def atom_dist(a, b):
    return math.sqrt((a.xx-b.xx)**2 + (a.xy-b.xy)**2 + (a.xz-b.xz)**2)

print("=" * 60)
print("exp-F Equilibrated Structure Analysis")
print("=" * 60)

print("\nLoading structures...")
sys_eq = pmd.load_file(PRMTOP, EQUIL)
sys_min = pmd.load_file(PRMTOP, MIN)
print(f"  Equilibrated: {len(sys_eq.atoms)} atoms, {len(sys_eq.residues)} residues")
print(f"  Minimized:    {len(sys_min.atoms)} atoms, {len(sys_min.residues)} residues")

# ============================================================================
# 1. Backbone RMSD (equil vs min) — NOTE: may be affected by PBC wrapping
# ============================================================================
print("\n--- 1. Backbone RMSD (equil vs min) ---")
print("  WARNING: NPT equilibration may cause PBC wrapping. RMSD values")
print("  below may be inflated if structures are in different images.")

# Simple RMSD without superposition (for diagnostic)
rec_ca_eq = [(a.xx, a.xy, a.xz) for r in sys_eq.residues[:1113] for a in r.atoms if a.name == 'CA']
rec_ca_min = [(a.xx, a.xy, a.xz) for r in sys_min.residues[:1113] for a in r.atoms if a.name == 'CA']
if len(rec_ca_eq) == len(rec_ca_min):
    d2 = sum((e[0]-m[0])**2 + (e[1]-m[1])**2 + (e[2]-m[2])**2 for e, m in zip(rec_ca_eq, rec_ca_min))
    print(f"  Receptor CA raw RMSD: {math.sqrt(d2/len(rec_ca_eq)):.3f} Å")

pep_ca_eq = [(a.xx, a.xy, a.xz) for r in sys_eq.residues[1113:1139] for a in r.atoms if a.name == 'CA']
pep_ca_min = [(a.xx, a.xy, a.xz) for r in sys_min.residues[1113:1139] for a in r.atoms if a.name == 'CA']
if len(pep_ca_eq) == len(pep_ca_min):
    d2 = sum((e[0]-m[0])**2 + (e[1]-m[1])**2 + (e[2]-m[2])**2 for e, m in zip(pep_ca_eq, pep_ca_min))
    print(f"  Peptide CA raw RMSD:  {math.sqrt(d2/len(pep_ca_eq)):.3f} Å")

# ============================================================================
# 2. LNK geometry
# ============================================================================
print("\n--- 2. LNK (Linker + C18) Geometry ---")
lnk = sys_eq.residues[1139]
lys = sys_eq.residues[1129]  # LYS P117

def get_atom(res, name):
    for a in res.atoms:
        if a.name == name:
            return a
    return None

nz = get_atom(lys, 'NZ')
c = get_atom(lnk, 'C')       # carbonyl C (attached to LYS NZ)
c13 = get_atom(lnk, 'C13')   # first CH2 after carbonyl
c55 = get_atom(lnk, 'C55')   # C18 diacid central C
o56 = get_atom(lnk, 'O56')   # carboxyl O1
o57 = get_atom(lnk, 'O57')   # carboxyl O2
c58 = get_atom(lnk, 'C58')   # gamma-Glu side chain CH2
o59 = get_atom(lnk, 'O59')   # gamma-Glu carboxyl O
n60 = get_atom(lnk, 'N60')   # gamma-Glu amide N
c61 = get_atom(lnk, 'C61')   # gamma-Glu CH3 (methyl end)

print(f"  LNK residue: {len(lnk.atoms)} atoms")
print(f"  LYS.NZ → LNK.C (amide bond):     {atom_dist(nz, c):.3f} Å")
print(f"  LNK.C → C13 (first CH2):         {atom_dist(c, c13):.3f} Å")
print(f"  LNK.C → C55 (C18 center):        {atom_dist(c, c55):.3f} Å")
print(f"  LNK.C → C61 (γGlu CH3):          {atom_dist(c, c61):.3f} Å")
print(f"  C55 → O56 (diacid):              {atom_dist(c55, o56):.3f} Å")
print(f"  C55 → O57 (diacid):              {atom_dist(c55, o57):.3f} Å")
print(f"  C58 → O59 (γGlu COOH):           {atom_dist(c58, o59):.3f} Å")
print(f"  C58 → N60 (γGlu amide N):        {atom_dist(c58, n60):.3f} Å")

# End-to-end: LYS NZ to farthest LNK heavy atom
heavy_atoms = [a for a in lnk.atoms if a.element != 1]
farthest = max(heavy_atoms, key=lambda a: atom_dist(nz, a))
print(f"  LYS.NZ → farthest heavy ({farthest.name}): {atom_dist(nz, farthest):.3f} Å")

# ============================================================================
# 3. Peptide-Receptor interface
# ============================================================================
print("\n--- 3. Peptide-Receptor Interface ---")

pep_atoms = [a for r in sys_eq.residues[1113:1139] for a in r.atoms]
rec_atoms = [a for r in sys_eq.residues[0:1113] for a in r.atoms]

pep_com = np.array([sum(a.xx for a in pep_atoms)/len(pep_atoms),
                    sum(a.xy for a in pep_atoms)/len(pep_atoms),
                    sum(a.xz for a in pep_atoms)/len(pep_atoms)])
rec_com = np.array([sum(a.xx for a in rec_atoms)/len(rec_atoms),
                    sum(a.xy for a in rec_atoms)/len(rec_atoms),
                    sum(a.xz for a in rec_atoms)/len(rec_atoms)])

print(f"  Peptide COM:  ({pep_com[0]:.2f}, {pep_com[1]:.2f}, {pep_com[2]:.2f})")
print(f"  Receptor COM: ({rec_com[0]:.2f}, {rec_com[1]:.2f}, {rec_com[2]:.2f})")
print(f"  COM separation: {np.linalg.norm(pep_com - rec_com):.2f} Å")

# Closest receptor contacts per peptide residue CA
print("\n  Peptide CA → closest receptor atom:")
pep_residues = sys_eq.residues[1113:1139]
for i, r in enumerate(pep_residues):
    ca = [a for a in r.atoms if a.name == 'CA']
    if not ca:
        continue
    ca = ca[0]
    closest = min(rec_atoms, key=lambda a: atom_dist(ca, a))
    print(f"    {r.name}{i+1:2d} → {closest.residue.name:3s} {closest.name:4s} (res {closest.residue.number}): {atom_dist(ca, closest):.2f} Å")

# ============================================================================
# 4. LNK tail extension and receptor proximity
# ============================================================================
print("\n--- 4. LNK Tail vs Receptor Proximity ---")

o56_closest = min(rec_atoms, key=lambda a: atom_dist(o56, a))
o57_closest = min(rec_atoms, key=lambda a: atom_dist(o57, a))
c61_closest = min(rec_atoms, key=lambda a: atom_dist(c61, a))

print(f"  LNK.O56 → nearest receptor: {o56_closest.residue.name} {o56_closest.name} (res {o56_closest.residue.number}): {atom_dist(o56, o56_closest):.2f} Å")
print(f"  LNK.O57 → nearest receptor: {o57_closest.residue.name} {o57_closest.name} (res {o57_closest.residue.number}): {atom_dist(o57, o57_closest):.2f} Å")
print(f"  LNK.C61 (γGlu CH3) → nearest: {c61_closest.residue.name} {c61_closest.name} (res {c61_closest.residue.number}): {atom_dist(c61, c61_closest):.2f} Å")

# Distance from C18 diacid to peptide COM
c55_pep = np.linalg.norm(np.array([c55.xx, c55.xy, c55.xz]) - pep_com)
print(f"  C55 (C18 center) → peptide COM: {c55_pep:.2f} Å")

# ============================================================================
# 5. LNK conformation: is it extended or folded?
# ============================================================================
print("\n--- 5. LNK Conformation Analysis ---")

# Distance from amide C to key points along the chain
c_to_c55 = atom_dist(c, c55)
c_to_c61 = atom_dist(c, c61)
c_to_farthest = atom_dist(c, farthest)

# For a fully extended 18-carbon chain + OEGs + gamma-Glu:
# C-C bond ~1.5 Å, C-C-C angle ~109.5°, each carbon adds ~1.25 Å projection
# 18C chain ≈ 22-25 Å extended
# 2×OEG (each ~7-8 Å) ≈ 14-16 Å
# gamma-Glu ≈ 5-6 Å
# Total extended ≈ 40-50 Å

print(f"  C (amide) → C55 (C18 center):    {c_to_c55:.1f} Å")
print(f"  C (amide) → C61 (γGlu end):      {c_to_c61:.1f} Å")
print(f"  C (amide) → farthest atom:       {c_to_farthest:.1f} Å")

if c_to_farthest > 35:
    print(f"  → Linker appears EXTENDED (>{c_to_farthest:.0f} Å end-to-end)")
elif c_to_farthest > 20:
    print(f"  → Linker appears PARTIALLY EXTENDED ({c_to_farthest:.0f} Å)")
else:
    print(f"  → Linker appears FOLDED/COMPACT ({c_to_farthest:.0f} Å)")

print("\n" + "=" * 60)
