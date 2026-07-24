#!/usr/bin/env python3
"""Pilot/production MD for exp-G: HSA + semaglutide + GLP-1R ECD ternary complex.

Differences from exp-D/md/run_md.py (physics unchanged otherwise):
- barostat added BEFORE Simulation() (same fix)
- weak CA restraints on HSA + ECD backbone ONLY (residues 1..682; peptide
  and LNK free), strength as a global parameter: 5.0 during heat/eq, 1.0
  in production
- anchor restraint on the FA3 salt-bridge oxygens (O56/O57/O38), k=5.0
  during heat/eq, released to 0.0 at production start so the competition
  is observed restraint-free
"""
import sys, os, argparse, time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_G = f"{REPO}/exps/exp-G"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gpu", type=str, default="3")
    p.add_argument("--nsteps", type=int, default=10_000_000)  # 20 ns pilot
    p.add_argument("--out", type=str, default="pilot")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()

def run(nsteps, gpu, out, seed):
    md_dir = f"{EXP_G}/md/{out}"
    os.makedirs(md_dir, exist_ok=True)
    prmtop = f"{EXP_G}/tleap/complex.prmtop"
    inpcrd = f"{EXP_G}/tleap/complex.inpcrd"

    print(f"Loading {prmtop}", flush=True)
    amber = pmd.load_file(prmtop, inpcrd)
    print(f"  Atoms: {len(amber.atoms)}, Residues: {len(amber.residues)}")

    system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometers,
                                constraints=app.HBonds, rigidWater=True)

    # CA restraints: HSA (pos 1-582) + ECD (583-682) backbone only, global kca
    # NOTE: parmed residue.number is 0-based on prmtop load -> HSA+ECD = 0..681
    ca_force = mm.CustomExternalForce("kca*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    ca_force.addGlobalParameter("kca", 5.0)
    ca_force.addPerParticleParameter("x0"); ca_force.addPerParticleParameter("y0")
    ca_force.addPerParticleParameter("z0")
    nca = 0
    for atom in amber.atoms:
        if atom.residue.number <= 681 and atom.name == "CA":
            xyz = amber.positions[atom.idx]
            ca_force.addParticle(atom.idx, [xyz[0], xyz[1], xyz[2]])
            nca += 1
    system.addForce(ca_force)
    print(f"CA restraints (HSA+ECD): {nca}")

    # anchor restraints on FA3 salt-bridge oxygens, global kanc
    anc_force = mm.CustomExternalForce("kanc*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    anc_force.addGlobalParameter("kanc", 5.0)
    anc_force.addPerParticleParameter("x0"); anc_force.addPerParticleParameter("y0")
    anc_force.addPerParticleParameter("z0")
    nanc = 0
    for atom in amber.atoms:
        if atom.residue.name == "LNK" and atom.name in ("O56", "O57", "O38"):
            xyz = amber.positions[atom.idx]
            anc_force.addParticle(atom.idx, [xyz[0], xyz[1], xyz[2]])
            nanc += 1
    system.addForce(anc_force)
    print(f"Anchor restraints (O56/O57/O38): {nanc}")

    # Barostat BEFORE Simulation()
    system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 310 * unit.kelvin))

    integrator = mm.LangevinIntegrator(310 * unit.kelvin, 1.0 / unit.picoseconds,
                                       2.0 * unit.femtoseconds)
    integrator.setRandomNumberSeed(seed)
    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(amber.topology, system, integrator, platform,
                                {"CudaDeviceIndex": gpu, "CudaPrecision": "mixed"})
    print("CUDA context ready", flush=True)

    simulation.reporters.append(app.DCDReporter(f"{md_dir}/traj.dcd", 50000))
    simulation.reporters.append(app.StateDataReporter(f"{md_dir}/log.txt", 10000,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True, speed=True))
    simulation.reporters.append(app.CheckpointReporter(f"{md_dir}/checkpoint.chk", 500000))

    simulation.context.setPositions(amber.positions)
    print("Minimizing...", flush=True)
    simulation.minimizeEnergy(maxIterations=10000)
    state = simulation.context.getState(getEnergy=True)
    print(f"  PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol", flush=True)

    print("Heating 0->100 K (NVT, 50 ps)"); integrator.setTemperature(100 * unit.kelvin); simulation.step(25000)
    print("Heating 100->310 K (NPT, 100 ps)")
    for i in range(5):
        integrator.setTemperature((100 + (i + 1) * 42) * unit.kelvin); simulation.step(10000)
    print("NPT eq (200 ps, 310 K, kca=5.0, kanc=5.0)")
    integrator.setTemperature(310 * unit.kelvin); simulation.step(100000)

    # graded release: CA 5.0->1.0, anchors fully off for production
    simulation.context.setParameter("kca", 1.0)
    simulation.context.setParameter("kanc", 0.0)
    print("Released: kca=1.0, kanc=0.0 -> production", flush=True)

    sim_start = time.time(); steps_done = 0
    while steps_done < nsteps:
        chunk = min(500000, nsteps - steps_done)
        simulation.step(chunk); steps_done += chunk
        elapsed = time.time() - sim_start
        ns_done = steps_done * 2e-6; ns_day = ns_done / (elapsed / 86400)
        remaining = (nsteps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] {ns_done:.0f}/{nsteps*2e-6:.0f} ns "
              f"({ns_day:.0f} ns/d, ~{remaining:.0f}h)", flush=True)
    print(f"Done! {((time.time() - sim_start) / 3600):.1f}h")

if __name__ == "__main__":
    args = parse_args()
    run(args.nsteps, args.gpu, args.out, args.seed)
