#!/usr/bin/env python3
"""
200 ns production MD for DPP-4 + GLP-1 peptide complexes (WT and Aib8).

Restraint strategy:
  - DPP-4 backbone CA: weak harmonic (10 kJ/mol/nm²), excluding active-site residues
  - Peptide N-term (residues 1-6 in chain P, = GLP-1 7-12): moderate backbone (N/CA/C)
    restraint (100 kJ/mol/nm²) to prevent long-peptide drift (CJC-1295 Trap 5)
  - Catalytic triad (S630, D708, H740), S1 pocket residues, and peptide scissile
    bond region (Aib8/Ala8, Glu9): NO restraints — free to show WT vs Aib8 differences

Protocol:
  - 0-100 K heating (NVT, 50 ps)
  - 100-310 K heating (NPT, 100 ps)
  - NPT equilibration (200 ps, 310 K, 1 bar)
  - NVT production (200 ns, 310 K)

Usage:
  python run_md.py --system wt   [--restart checkpoint]
  python run_md.py --system aib8 [--restart checkpoint]
"""
import sys, os, argparse, time, signal
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_A = f"{REPO}/exps/exp-A"

# Active-site residues to EXCLUDE from DPP-4 CA restraint
# Catalytic triad: S630, D708, H740
# S1 pocket: Y547, W629, Y631, V656, Y662, Y666, N710 (near S1 pocket)
# Also exclude residues within 5A of catalytic Ser630
DPP4_FREE_RESIDUES = {630, 708, 740, 547, 629, 631, 656, 662, 666, 710, 711, 631, 632, 633, 634, 635, 547, 548}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--system", required=True, choices=["wt", "aib8"])
    p.add_argument("--restart", default=None, help="Path to checkpoint XML")
    p.add_argument("--nsteps", type=int, default=None, help="Override production steps")
    return p.parse_args()

def add_restraints(system, topology, positions):
    """Add harmonic position restraints to the system."""
    # DPP-4 CA restraint (weak, 10 kJ/mol/nm²)
    dpp4_ca_force = mm.CustomExternalForce("10.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    dpp4_ca_force.addPerParticleParameter("x0")
    dpp4_ca_force.addPerParticleParameter("y0")
    dpp4_ca_force.addPerParticleParameter("z0")

    # Peptide N-term backbone restraint (moderate, 100 kJ/mol/nm²)
    pep_bb_force = mm.CustomExternalForce("100.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    pep_bb_force.addPerParticleParameter("x0")
    pep_bb_force.addPerParticleParameter("y0")
    pep_bb_force.addPerParticleParameter("z0")

    # Identify atoms
    dpp4_restrained = 0
    pep_restrained = 0
    for chain in topology.chains():
        for res in chain.residues():
            # Peptide chain P: restrain N-term residues 1-6 (GLP-1 7-12) backbone N/CA/C
            if chain.id == 'P':
                if int(res.id) <= 6:
                    for atom in res.atoms():
                        if atom.name in ('N', 'CA', 'C'):
                            xyz = positions[atom.index]
                            pep_bb_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                            pep_restrained += 1
                continue

            # DPP-4 chain: CA restraint excluding active site
            if res.name in ('WAT', 'Na+', 'Cl-', 'HOH', 'SOL'):
                continue
            try:
                resid = int(res.id)
            except ValueError:
                continue
            if resid in DPP4_FREE_RESIDUES:
                continue
            for atom in res.atoms():
                if atom.name == 'CA':
                    xyz = positions[atom.index]
                    dpp4_ca_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                    dpp4_restrained += 1

    idx = system.addForce(dpp4_ca_force)
    idx = system.addForce(pep_bb_force)
    print(f"Restraints: {dpp4_restrained} DPP-4 CA + {pep_restrained} peptide N-term BB")
    return system

def run(mode, restart=None, nsteps=None):
    """Run MD for a given system."""
    md_dir = f"{EXP_A}/md/{mode}"
    system_xml = f"{md_dir}/{mode}_system.xml"
    pdb_path = f"{md_dir}/{mode}_minimized.pdb"

    if not os.path.exists(system_xml):
        # Rebuild system from prmtop
        prmtop = f"{md_dir}/{mode}_modified.prmtop" if mode == "aib8" else f"{EXP_A}/tleap/{mode}.prmtop"
        inpcrd = f"{md_dir}/{mode}_modified.inpcrd" if mode == "aib8" else f"{EXP_A}/tleap/{mode}.inpcrd"
        if not os.path.exists(prmtop):
            prmtop = f"{EXP_A}/tleap/{mode}.prmtop"
            inpcrd = f"{EXP_A}/tleap/{mode}.inpcrd"
        print(f"Rebuilding system from {prmtop}")
        amber = pmd.load_file(prmtop, inpcrd)
        system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                                     constraints=app.HBonds, rigidWater=True)
        with open(system_xml, "w") as f:
            f.write(mm.XmlSerializer.serialize(system))
        # Save PDB for topology
        amber.save(f"{md_dir}/{mode}_amber.pdb", overwrite=True)
        pdb_path = f"{md_dir}/{mode}_amber.pdb"

    # Load system
    print(f"Loading system from {system_xml}")
    with open(system_xml) as f:
        system = mm.XmlSerializer.deserialize(f.read())

    # Load topology
    pdb = app.PDBFile(pdb_path)

    # Add restraints
    system = add_restraints(system, pdb.topology, pdb.positions)

    # Integrator: 2 fs, 310 K, Langevin
    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)

    # Platform: CUDA
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

    # Set up checkpointing and reporting
    report_interval = 10000  # every 20 ps
    checkpoint_interval = 500000  # every 1 ns
    dcd_interval = 50000  # every 100 ps

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
        simulation.context.setPositions(pdb.positions)

    # Determine production steps
    if nsteps is None:
        nsteps = 100_000_000  # 200 ns at 2 fs = 100M steps

    # Heating protocol (skip if restarting)
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

    # Track total steps
    sim_start = time.time()
    steps_done = 0
    report_next = 0

    while steps_done < nsteps:
        chunk = min(500000, nsteps - steps_done)  # 1 ns chunks
        simulation.step(chunk)
        steps_done += chunk

        # Progress report
        elapsed = time.time() - sim_start
        ns_done = steps_done * 2e-6
        ns_day = ns_done / (elapsed / 86400)
        remaining = (nsteps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] {ns_done:.0f}/{nsteps*2e-6:.0f} ns "
              f"({ns_day:.0f} ns/day, ~{remaining:.0f}h remaining)")

        # Save checkpoint at end
        simulation.saveCheckpoint(f"{md_dir}/{mode}_final.chk")

    print(f"Done! Total time: {(time.time()-sim_start)/3600:.1f}h")

if __name__ == "__main__":
    args = parse_args()
    run(args.system, args.restart, args.nsteps)
