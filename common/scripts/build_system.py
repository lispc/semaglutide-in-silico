#!/usr/bin/env python3
"""
Build, solvate, and minimize DPP-4 + GLP-1 peptide complex using OpenMM + ff14SB.

Usage:
  python build_system.py --peptide wt
  python build_system.py --peptide aib8
"""
import sys, os, argparse, copy
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_A = f"{REPO}/exps/exp-A"

def build_and_minimize(peptide_label):
    pdb_path = f"{EXP_A}/structures/DPP4_GLP1_{'WT' if peptide_label == 'wt' else 'Aib8'}_clean.pdb"
    out_prefix = f"{EXP_A}/md/{peptide_label}/{peptide_label}"
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    print(f"=== Building {peptide_label} system ===")

    pdb = app.PDBFile(pdb_path)

    forcefield = app.ForceField(
        'amber14/protein.ff14SB.xml',
        'amber14/tip3p.xml',
        f'{REPO}/common/params/aib_residue.xml'
    )
    print("Force field: ff14SB + TIP3P + Aib custom")

    modeller = app.Modeller(pdb.topology, pdb.positions)

    # Try adding hydrogens; if it fails due to terminal residue template
    # mismatches, delete that residue from the topology and retry.
    # This handles PDB structures with non-standard termini or incomplete chains.
    max_retries = 10
    for attempt in range(max_retries):
        try:
            modeller.addHydrogens(forcefield, pH=7.0)
            break
        except ValueError as e:
            err = str(e)
            # Parse the residue index from the error: "residue %d (XXX)"
            import re
            m = re.search(r'residue (\d+)', err)
            if not m:
                raise
            bad_idx = int(m.group(1))
            # Find this residue in topology
            to_delete = []
            for chain in modeller.topology.chains():
                for res in chain.residues():
                    if res.index == bad_idx:
                        to_delete.append(res)
                        break
            if not to_delete:
                # Try deleting by name instead - get the residue with this index
                for chain in modeller.topology.chains():
                    residues = list(chain.residues())
                    # Try the last residue of the chain with the most residues
                    pass
                raise
            modeller.delete(to_delete)
            print(f"  Deleted residue idx {bad_idx} ({to_delete[0].name}) — template mismatch, retrying...")

    modeller.addSolvent(forcefield, padding=1.2*unit.nanometers,
                        model='tip3p', ionicStrength=0.1*unit.molar)
    na = modeller.topology.getNumAtoms()
    print(f"System: {na} atoms")

    with open(f"{out_prefix}_solvated.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
        ignoreExternalBonds=True
    )

    # Minimize
    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    pe = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Initial PE: {pe}")
    simulation.minimizeEnergy(maxIterations=5000)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Minimized PE: {state.getPotentialEnergy()}")

    with open(f"{out_prefix}_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    with open(f"{out_prefix}_system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    print(f"Done: {out_prefix}_*\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--peptide", required=True, choices=["wt", "aib8"])
    args = parser.parse_args()
    build_and_minimize(args.peptide)
