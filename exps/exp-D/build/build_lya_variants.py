#!/usr/bin/env python3
"""Build 5 LYA linker variants for exp-D (Lau 2015 Table 3).
Pipeline: RDKit 3D → mol2 → antechamber AM1-BCC → positioned mol2
"""
import os, sys, subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

OUT_DIR = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"

# C18 diacid tail: -(CH2)16-COO⁻ (16 CH2 + 1 COO⁻ = 17 carbons here,
# plus the amide/GGlu carbonyl = 18C total for the diacid)
TAIL = "CCCCCCCCCCCCCCCCC(=O)[O-]"

# Lys side chain (without the linker part): CCCC = -(CH2)4-
LYS_SIDE = "CCCC"

# Variant definitions: (name, linker_smiles_fragment, description)
# Pattern: ...LYS_SIDE-NC(=O)-{linker_frag}-{TAIL}
# The NC(=O) is hardcoded as it represents the Lys-NZ amide bond
# For no_linker: linker_frag is empty (TAIL directly follows NC(=O))
VARIANTS = [
    ("no_linker",  "",                                "Lys26直达C18二酸 (Cmpd 19)"),
    ("gglu",       "CCC(N)C(=O)",                     "γGlu linker (Cmpd 20)"),
    ("gglu_1oeg",  "CCC(N)C(=O)NCCOCCOCC(=O)",        "γGlu-1×OEG (Cmpd 21)"),
    ("gglu_2oeg",  "CCC(N)C(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)",
                                                      "γGlu-2×OEG (Semaglutide)"),
    ("gglu_3oeg",  "CCC(N)C(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)",
                                                      "γGlu-3×OEG (Cmpd 24)"),
]

def build_lya(variant_name, linker_frag, description, net_charge=0):
    """Build one LYA variant."""
    # Full SMILES: ACE - Lys(backbone) - Lys side chain - linker - C18 tail - NME
    smiles = f"CC(=O)NC({LYS_SIDE}NC(=O){linker_frag}{TAIL})C(=O)NC"
    os.makedirs(OUT_DIR, exist_ok=True)
    os.chdir(OUT_DIR)

    print(f"\n{'='*60}")
    print(f"Building: {variant_name} ({description})")
    print(f"SMILES: {smiles[:80]}...")

    # Step 1: RDKit 3D
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("  ERROR: Invalid SMILES")
        return False
    mol = Chem.AddHs(mol)
    print(f"  Atoms: {mol.GetNumAtoms()} (before embedding)")

    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        print(f"  ETKDGv3 returned {result}, trying random coords...")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    # Step 2: Write raw mol2 WITHOUT atom types (empty column)
    # This lets antechamber assign GAFF2 types from scratch
    raw_mol2 = f"{variant_name}_raw.mol2"

    with open(raw_mol2, 'w') as f:
        n = mol.GetNumAtoms()
        f.write("@<TRIPOS>MOLECULE\nLYA\n")
        f.write(f" {n} 0 0 0 0\nSMALL\nGAFF2\n{description}\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i in range(n):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            name = f"{atom.GetSymbol()}{i+1}"
            # Empty atom type column → antechamber assigns GAFF2 type
            f.write(f"{i+1:6d} {name:6s} {pos.x:10.4f} {pos.y:10.4f} {pos.z:10.4f} {name[:2]:4s} 1 LYA 0.0000\n")
        f.write("@<TRIPOS>BOND\n")
        bid = 0
        for bond in mol.GetBonds():
            bid += 1
            f.write(f"{bid:6d} {bond.GetBeginAtomIdx()+1:6d} {bond.GetEndAtomIdx()+1:6d} {int(bond.GetBondTypeAsDouble())}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 LYA 1 TEMP 0 **** **** 0 ROOT\n")

    print(f"  {raw_mol2} written ({n} atoms)")

    # Step 3: antechamber AM1-BCC
    bcc_mol2 = f"{variant_name}_bcc.mol2"
    print(f"  Running AM1-BCC (nc={net_charge})...")
    cmd = (f"rm -f ANTECHAMBER* ATOMTYPE* sqm.* 2>/dev/null && "
           f"antechamber -i {raw_mol2} -fi mol2 -o {bcc_mol2} -fo mol2 "
           f"-c bcc -at gaff2 -nc {net_charge} -rn LYA -pf yes 2>&1")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                            timeout=3600, executable='/bin/bash')
    out = result.stdout + result.stderr
    if "Fatal" in out or "Error" in out or result.returncode != 0:
        print(f"  FAILED:\n{out[-500:]}")
        return False

    # Check charges
    charges = []
    if os.path.exists(bcc_mol2):
        with open(bcc_mol2) as f:
            in_atom = False
            for line in f:
                if "@<TRIPOS>ATOM" in line: in_atom = True; continue
                if "@<TRIPOS>BOND" in line: break
                if in_atom and line.strip():
                    parts = line.split()
                    if len(parts) >= 9: charges.append(float(parts[-1]))
        print(f"  AM1-BCC OK: {len(charges)} atoms, total={sum(charges):.4f}, range=[{min(charges):.3f}, {max(charges):.3f}]")
    else:
        print(f"  WARNING: {bcc_mol2} not found")
        return False

    return True

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    success = 0
    for name, linker, desc in VARIANTS:
        # All variants: [O-] terminal → net formal charge = -1
        # nc=-1 corrects e⁻ count to even for SQM convergence
        nc = -1
        if build_lya(name, linker, desc, nc):
            success += 1
    print(f"\n{'='*60}")
    print(f"Done: {success}/{len(VARIANTS)} built successfully")
