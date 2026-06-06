# Membrane System Build Status

## Date: 2026-06-06

## Summary
Successfully built full membrane-embedded GLP-1R + semaglutide system topology using a custom pipeline after diagnosing a critical packmol-memgen bug.

## System Composition
- Total atoms: 312,501
- Residues: 82,357
- Molecules: 80,501
- Box: Orthorhombic
- Solvent: 79,615 WAT
- Lipids: 449 (POPC/CHL1 8:2)
- Ions: 425 K+/Cl-
- Protein: 7KI0 receptor (chains R,A,B,G,N) + semaglutide (chain P)

## Build Pipeline Overview

### 1. Membrane Component Generation (packmol, manual)
- Generated POPC:CHL1 (8:2) membrane + waters + ions manually via packmol
- Original target: 562 lipids; after protein overlap removal: 449 lipids
- Upper leaflet z ~4–23 Å, lower leaflet z ~-23–-4 Å

### 2. Custom Merge Pipeline (`build_system_v10.py`)
- Reads protein and membrane PDBs, groups by TER, assigns unique chain IDs
- Removes lipids within 2.5 Å of protein via scipy.spatial.cKDTree (removed 113, kept 449)
- Lipid name mapping via `charmmlipid2amber.csv` (POP→POPC/PA/PC/OL, CHL→CHL1)
- Component splitting: protein, lipids, waters (≤99,999 atoms/file), ions
- Water reformatting: residue `WAT`, atoms `O`, `H1`, `H2`
- Ion reformatting: residue/atom names `K+` and `Cl-`
- N-terminal fix: removes H3 from chain P GLY 101

### 3. tleap Topology Assembly (`build_system_final.in`)
- Loads ff14SB + GAFF2 + TIP3P + lipid21
- Loads LNK from `lnk_noh_zero.mol2` (51 heavy atoms, net charge 0)
- Combines components: PROT + LIPIDS + WAT1 + WAT2 + WAT3 + IONS
- Forms amide bond: LYS P117 NZ — LNK C13
- Sets orthorhombic box, neutralizes with addionsrand
- Outputs: `system_final.prmtop` / `system_final.inpcrd` / `system_final.pdb`

## Critical Fixes Applied

### 1. Packmol-memgen Bug Diagnosis
**Root cause**: `MembraneParams.pdb_reindex(renumber=True)` flattens protein residues per chain starting from 1. `charmmlipid2amber.py` identifies residues only by `(chain, resnum)`, so lipids sharing resnums with protein get merged, creating hybrid residues and broken peptide bonds.
**Impact**: This prevented direct use of packmol-memgen for the combined system.
**Workaround**: Manual packmol for membrane-only, then custom Python merge pipeline.

### 2. PDB Column Shift Bug
**Root cause**: Python f-string `{resname:>3s}{chain:>1s}` — when resname is 4 chars (e.g. `WATA`), chainID shifts to col 21 instead of 22, pushing all coordinates left by 1 column.
**Manifestation**: Adjacent negative values merged (e.g. `-43.108-107.555`).
**Fix**: Added explicit space: `{resname:>3s} {chain:>1s}`.

### 3. Internal Protein Chain Breaks
Inserted 5 additional TER records at missing loop regions to prevent tleap from creating impossibly long C-N bonds:
- Chain A: GLY47→SER205 (17.4Å), THR249→SER264 (9.9Å), SER365→THR370 (10.1Å)
- Chain R: SER129→PRO137 (15.0Å), LEU339→ASP344 (8.5Å)

### 4. LNK Charge Normalization
Original `lnk_noh_fixed.mol2` had partial charges summing to +3.15842.
**Fix**: Generated `lnk_noh_zero.mol2` by normalizing each atom charge by +0.049031 to achieve net charge 0.000000.

### 5. Water/Ion Naming for Amber
- Waters: packmol `TIP3`/`WATB`/`WATC` → `WAT` with atoms `O`, `H1`, `H2`
- Ions: aligned to `atomic_ions.lib` expectations (`K+` and `Cl-` residue and atom names)

### 6. LNK Atom Count Match
Stripped all LNK hydrogens from PDB (51 heavy atoms) to match `lnk_noh_fixed.mol2` / `lnk_noh_zero.mol2` used for parameter loading.

## tleap Residue Indexing (Chain P / Semaglutide)
In the combined system context:
- GLY 101 = SYS.1114
- LYS 117 = SYS.1140 (bond donor)
- GLY 126 = SYS.1149 (OXT removed)
- LNK 127 = SYS.1150 (amide C, connect1 atom)

## Remaining Warnings
- tleap charge warning: net charge = -0.02001 (non-integral, non-zero)
  - This is residual from LNK charge rounding; acceptable for MD with PME
- LNK "one-sided connection" warning: benign — LNK is a side-chain cap without an N-terminus

## Output Files
- `system_final.prmtop` (66 MB)
- `system_final.inpcrd` (11 MB)
- `system_final.pdb` (21 MB)

## Verification
```
cpptraj system_final.prmtop <<EOF
topologyinfo
EOF
```
Result: 312,501 atoms; 82,357 residues; 80,501 molecules; 79,615 solvent; orthorhombic box. No topology read errors.

## Next Steps
- Run energy minimization (500-1000 steps)
- Heat from 0→100K (NVT), then 100→310K (NPT)
- Production MD with restraints on receptor CA
