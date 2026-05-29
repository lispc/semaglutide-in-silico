#!/usr/bin/env python3
"""
Production MD for HSA + fatty acid at FA3 site. Supports replicas.

Usage:
  python run_md.py --system c18_monoacid --replica 1 --gpu 0
  python run_md.py --system c18_monoacid --replica 2 --gpu 1
"""
import sys, os, argparse, time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_C = f"{REPO}/exps/exp-C"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--system", required=True)
    p.add_argument("--replica", type=int, default=1)
    p.add_argument("--restart", default=None)
    p.add_argument("--nsteps", type=int, default=50_000_000)  # 100 ns default
    p.add_argument("--gpu", type=str, default="0")
    return p.parse_args()

def add_restraints(system, amber):
    """Add weak CA restraint on HSA backbone."""
    ca_force = mm.CustomExternalForce("5.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    ca_force.addPerParticleParameter("x0"); ca_force.addPerParticleParameter("y0"); ca_force.addPerParticleParameter("z0")
    count = 0
    for atom in amber.atoms:
        if atom.residue.name not in ('WAT','HOH','SOL','Na+','Cl-','FAH') and atom.name == 'CA':
            xyz = amber.positions[atom.idx]
            ca_force.addParticle(atom.idx, [xyz[0], xyz[1], xyz[2]])
            count += 1
    system.addForce(ca_force)
    print(f"Restraints: {count} HSA CA")
    return system

def run(system_name, replica=1, restart=None, nsteps=50_000_000, gpu="0"):
    md_dir = f"{EXP_C}/md/{system_name}/rep{replica}"
    os.makedirs(md_dir, exist_ok=True)

    prmtop = f"{EXP_C}/tleap/{system_name}.prmtop"
    inpcrd = f"{EXP_C}/tleap/{system_name}.inpcrd"

    print(f"Loading {prmtop}")
    amber = pmd.load_file(prmtop, inpcrd)
    system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds, rigidWater=True)
    system = add_restraints(system, amber)

    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    integrator.setRandomNumberSeed(replica * 42)
    platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(amber.topology, system, integrator, platform,
                                {'CudaDeviceIndex': gpu, 'CudaPrecision': 'mixed'})

    simulation.reporters.append(app.DCDReporter(f"{md_dir}/{system_name}_traj.dcd", 50000))
    simulation.reporters.append(app.StateDataReporter(f"{md_dir}/{system_name}_log.txt", 10000,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True))
    simulation.reporters.append(app.CheckpointReporter(f"{md_dir}/{system_name}_checkpoint.chk", 500000))

    if restart:
        print(f"Restarting from {restart}")
        with open(restart, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
    else:
        simulation.context.setPositions(amber.positions)
        # Minimization to resolve clashes
        print("Minimizing (10k steps)...")
        simulation.minimizeEnergy(maxIterations=10000)
        state = simulation.context.getState(getEnergy=True)
        pe_min = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        print(f"  PE after minimization: {pe_min:.0f} kJ/mol")

        # Graduated heating with strong→weak position restraints
        for k_val, label, temp, steps in [
            (200.0, "0→10K",   10,  5000),
            (200.0, "10→100K", 100, 25000),
            (100.0, "100→200K",200, 25000),
            (50.0,  "200→310K",310, 25000),
            (20.0,  "310K relax1",310,25000),
            (10.0,  "310K relax2",310,25000),
        ]:
            rf = mm.CustomExternalForce(f"{k_val} * periodicdistance(x, y, z, x0, y0, z0)^2")
            rf.addPerParticleParameter("x0"); rf.addPerParticleParameter("y0"); rf.addPerParticleParameter("z0")
            for atom in amber.atoms:
                if atom.residue.name not in ('WAT','HOH','SOL','Na+','Cl-'):
                    xyz = amber.positions[atom.idx]
                    rf.addParticle(atom.idx, [xyz[0], xyz[1], xyz[2]])
            rf_idx = system.addForce(rf)
            integrator.setTemperature(temp * unit.kelvin)
            print(f"  {label} ({k_val:.0f} kJ/mol/nm², {steps*2e-6*1000:.0f} ps)")
            simulation.step(steps)
            system.removeForce(rf_idx)

        print("NPT eq (100 ps, with barostat)")
        system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))
        integrator.setTemperature(310)
        simulation.step(50000)

    print(f"Production: {nsteps} steps ({nsteps*2e-6:.0f} ns)")
    sim_start = time.time(); steps_done = 0
    while steps_done < nsteps:
        chunk = min(500000, nsteps - steps_done)
        simulation.step(chunk); steps_done += chunk
        elapsed = time.time() - sim_start
        ns_done = steps_done * 2e-6; ns_day = ns_done / (elapsed / 86400)
        remaining = (nsteps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] {ns_done:.0f}/{nsteps*2e-6:.0f} ns ({ns_day:.0f} ns/d, ~{remaining:.0f}h)")

    simulation.saveCheckpoint(f"{md_dir}/{system_name}_final.chk")
    print(f"Done! {((time.time()-sim_start)/3600):.1f}h")

if __name__ == "__main__":
    args = parse_args()
    run(args.system, args.replica, args.restart, args.nsteps, args.gpu)
