#!/usr/bin/env python3
"""
Minimization + equilibration for exp-F: full semaglutide + GLP-1R receptor.

Phase 1 goal: validate system stability after solvation.

Restraint strategy:
  - Receptor backbone CA (residues 0-1112): weak harmonic (10 kJ/mol/nm²)
  - Peptide N-term (residues 1113-1118) backbone N/CA/C: moderate (100 kJ/mol/nm²)
  - LNK and peptide body: free

Protocol:
  - Minimization (500 steps L-BFGS)
  - 0→100 K heating (NVT, 50 ps)
  - 100→310 K heating (NPT, 100 ps)
  - NPT equilibration (100 ps, 310 K, 1 bar) — short for validation
  - Save final coordinates and checkpoint

Usage:
  python run_min_equil.py [--min-only]
"""
import sys, os, argparse, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/tleap/system.prmtop"
INPCRD = f"{EXP_F}/tleap/system.inpcrd"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--min-only", action="store_true", help="Only run minimization")
    p.add_argument("--restart", default=None, help="Path to checkpoint XML")
    return p.parse_args()

def add_restraints(system, topology, positions):
    """Add harmonic position restraints."""
    receptor_force = mm.CustomExternalForce("10.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    receptor_force.addPerParticleParameter("x0")
    receptor_force.addPerParticleParameter("y0")
    receptor_force.addPerParticleParameter("z0")

    pep_force = mm.CustomExternalForce("100.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    pep_force.addPerParticleParameter("x0")
    pep_force.addPerParticleParameter("y0")
    pep_force.addPerParticleParameter("z0")

    receptor_count = 0
    pep_count = 0
    res_idx = 0
    for chain in topology.chains():
        for res in chain.residues():
            if res.name in ('WAT', 'Na+', 'Cl-', 'HOH', 'SOL', 'LNK'):
                res_idx += 1
                continue
            # Peptide residues: 1113-1118 (N-term) backbone restraint
            if 1113 <= res_idx <= 1118:
                for atom in res.atoms():
                    if atom.name in ('N', 'CA', 'C'):
                        xyz = positions[atom.index]
                        pep_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                        pep_count += 1
                res_idx += 1
                continue
            # Receptor residues: 0-1112 CA restraint
            if res_idx <= 1112:
                for atom in res.atoms():
                    if atom.name == 'CA':
                        xyz = positions[atom.index]
                        receptor_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                        receptor_count += 1
            res_idx += 1

    system.addForce(receptor_force)
    system.addForce(pep_force)
    print(f"Restraints: {receptor_count} receptor CA + {pep_count} peptide N-term BB")
    return system

def run(min_only=False, restart=None):
    md_dir = f"{EXP_F}/md/rep1"
    os.makedirs(md_dir, exist_ok=True)

    print(f"Loading Amber system from {PRMTOP}")
    amber = pmd.load_file(PRMTOP, INPCRD)
    system = amber.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds,
        rigidWater=True,
    )

    topology = app.AmberPrmtopFile(PRMTOP).topology
    positions = amber.positions

    system = add_restraints(system, topology, positions)

    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

    system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))

    simulation = app.Simulation(topology, system, integrator, platform, properties)

    # Add reporter for progress monitoring
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, 5000,
        step=True, time=True, potentialEnergy=True, temperature=True,
        volume=True, speed=True
    ))

    if restart:
        print(f"Restarting from {restart}")
        with open(restart, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
    else:
        simulation.context.setPositions(positions)

    if not restart:
        print("Minimizing energy...")
        start = time.time()
        simulation.minimizeEnergy(maxIterations=500)
        print(f"Minimization done in {time.time()-start:.1f}s")

        state = simulation.context.getState(getPositions=True)
        with open(f"{md_dir}/minimized.pdb", "w") as f:
            app.PDBFile.writeFile(topology, state.getPositions(), f)

        if min_only:
            return

    if not restart:
        print("Heating 0→100 K (NVT, 50 ps)")
        integrator.setTemperature(100)
        simulation.step(25000)

        print("Heating 100→310 K (NPT, 100 ps)")
        for i in range(5):
            integrator.setTemperature(100 + (i+1)*42)
            simulation.step(20000)

        print("NPT equilibration (100 ps, 310 K)")
        integrator.setTemperature(310)
        simulation.step(50000)

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(f"{md_dir}/equilibrated.pdb", "w") as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f)
    simulation.saveCheckpoint(f"{md_dir}/equil.chk")
    print("Saved equilibrated.pdb and checkpoint")

    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"Final potential energy: {pe:.1f} kJ/mol")

if __name__ == "__main__":
    args = parse_args()
    run(min_only=args.min_only, restart=args.restart)
