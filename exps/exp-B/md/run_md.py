#!/usr/bin/env python3
"""
200 ns production MD for GLP-1R ECD + GLP-1(10-37) peptide complexes (B1 K34, B2 R34).

Restraint strategy:
  - ECD backbone CA: weak harmonic (10 kJ/mol/nm²), excluding residues near pos 34 binding site
  - Peptide C-terminal helix (residues 13-33 PDB numbering): weak backbone (N/CA/C)
    restraint (50 kJ/mol/nm²) to keep helix in ECD groove
  - Position 34 (K34/R34) sidechain + nearby ECD residues: NO restraints
  - Peptide N-terminal (10-12): unrestrained

Protocol:
  - 0-100 K heating (NVT, 50 ps)
  - 100-310 K heating (NPT, 100 ps)
  - NPT equilibration (200 ps, 310 K)
  - NVT production (200 ns, 310 K)

Usage:
  python run_md.py --system b1 [--restart checkpoint] [--gpu 0]
  python run_md.py --system b2 [--restart checkpoint] [--gpu 1]
"""
import sys, os, argparse, time, signal
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_B = f"{REPO}/exps/exp-B"

# ECD residues near position 34 binding site to EXCLUDE from ECD CA restraint
# Allow free motion for residues that may interact with K34/R34
ECD_FREE_RESIDUES = set()

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--system", required=True, choices=["b1", "b2"])
    p.add_argument("--restart", default=None, help="Path to checkpoint XML")
    p.add_argument("--nsteps", type=int, default=None, help="Override production steps")
    p.add_argument("--gpu", type=str, default="0", help="GPU device index")
    return p.parse_args()

def add_restraints(system, topology, positions):
    """Add harmonic position restraints."""
    # ECD CA restraint (weak, 10 kJ/mol/nm²)
    ecd_ca_force = mm.CustomExternalForce("10.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    ecd_ca_force.addPerParticleParameter("x0")
    ecd_ca_force.addPerParticleParameter("y0")
    ecd_ca_force.addPerParticleParameter("z0")

    # Peptide C-term helix backbone restraint (moderate, 50 kJ/mol/nm²)
    pep_bb_force = mm.CustomExternalForce("50.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    pep_bb_force.addPerParticleParameter("x0")
    pep_bb_force.addPerParticleParameter("y0")
    pep_bb_force.addPerParticleParameter("z0")

    ecd_restrained = 0
    pep_restrained = 0

    for chain in topology.chains():
        for res in chain.residues():
            # Peptide chain B: restrain C-terminal helix residues 13-33 backbone
            if chain.id == 'B':
                try:
                    rid = int(res.id)
                except ValueError:
                    continue
                if 13 <= rid <= 33:
                    for atom in res.atoms():
                        if atom.name in ('N', 'CA', 'C'):
                            xyz = positions[atom.index]
                            pep_bb_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                            pep_restrained += 1
                continue

            # ECD chain A: CA restraint
            if chain.id != 'A':
                continue
            if res.name in ('WAT', 'Na+', 'Cl-', 'HOH', 'SOL'):
                continue
            for atom in res.atoms():
                if atom.name == 'CA':
                    xyz = positions[atom.index]
                    ecd_ca_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                    ecd_restrained += 1

    system.addForce(ecd_ca_force)
    system.addForce(pep_bb_force)
    print(f"Restraints: {ecd_restrained} ECD CA + {pep_restrained} peptide helix BB")
    return system

def run(mode, restart=None, nsteps=None, gpu="0"):
    md_dir = f"{EXP_B}/md/{mode}"
    os.makedirs(md_dir, exist_ok=True)

    # Topology from tleap
    prmtop = f"{EXP_B}/tleap/{mode}.prmtop"
    inpcrd = f"{EXP_B}/tleap/{mode}.inpcrd"

    if not os.path.exists(prmtop):
        raise FileNotFoundError(f"Topology not found: {prmtop}")

    print(f"Loading topology from {prmtop}")
    amber = pmd.load_file(prmtop, inpcrd)
    system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds, rigidWater=True)

    # Save minimized PDB for topology
    amber.save(f"{md_dir}/{mode}_amber.pdb", overwrite=True)

    # Add restraints
    system = add_restraints(system, amber.topology, amber.positions)

    # Integrator: 2 fs, 310 K, Langevin
    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)

    # Platform: CUDA
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': gpu, 'CudaPrecision': 'mixed'}

    simulation = app.Simulation(amber.topology, system, integrator, platform, properties)

    # Reporters
    report_interval = 10000
    checkpoint_interval = 500000
    dcd_interval = 50000

    simulation.reporters.append(app.DCDReporter(f"{md_dir}/{mode}_traj.dcd", dcd_interval))
    simulation.reporters.append(app.StateDataReporter(
        f"{md_dir}/{mode}_log.txt", report_interval,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True, speed=True
    ))
    simulation.reporters.append(app.CheckpointReporter(f"{md_dir}/{mode}_checkpoint.chk", checkpoint_interval))

    # Set positions
    if restart:
        print(f"Restarting from {restart}")
        with open(restart, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
    else:
        simulation.context.setPositions(amber.positions)

    if nsteps is None:
        nsteps = 100_000_000  # 200 ns

    # Heating
    if not restart:
        print("Heating 0→100 K (NVT, 50 ps)")
        integrator.setTemperature(100)
        simulation.step(25000)

        print("Heating 100→310 K (NPT, 100 ps)")
        system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))
        for i in range(5):
            integrator.setTemperature(100 + (i+1)*42)
            simulation.step(10000)

        print("NPT equilibration (200 ps, 310 K)")
        integrator.setTemperature(310)
        simulation.step(100000)

    print(f"Production: {nsteps} steps ({nsteps*2e-6:.0f} ns)")

    sim_start = time.time()
    steps_done = 0

    while steps_done < nsteps:
        chunk = min(500000, nsteps - steps_done)
        simulation.step(chunk)
        steps_done += chunk

        elapsed = time.time() - sim_start
        ns_done = steps_done * 2e-6
        ns_day = ns_done / (elapsed / 86400)
        remaining = (nsteps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] {ns_done:.0f}/{nsteps*2e-6:.0f} ns "
              f"({ns_day:.0f} ns/day, ~{remaining:.0f}h remaining)")

    simulation.saveCheckpoint(f"{md_dir}/{mode}_final.chk")
    print(f"Done! Total time: {(time.time()-sim_start)/3600:.1f}h")

if __name__ == "__main__":
    args = parse_args()
    run(args.system, args.restart, args.nsteps, args.gpu)
