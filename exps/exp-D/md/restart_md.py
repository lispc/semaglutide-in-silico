#!/usr/bin/env python3
"""Restart MD from checkpoint. Usage: python restart_md.py --system no_linker_fixed --replica 1 --gpu 0"""
import sys, os, argparse, time, numpy as np
import openmm as mm, openmm.app as app, openmm.unit as unit, parmed as pmd

EXP_D = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D"

p = argparse.ArgumentParser()
p.add_argument("--system", required=True); p.add_argument("--replica", type=int, default=1)
p.add_argument("--gpu", type=str, default="0"); p.add_argument("--nsteps", type=int, default=50_000_000)
args = p.parse_args()

vname = args.system; rep = args.replica
md_dir = f"{EXP_D}/md/{vname}/rep{rep}"
prmtop = f"{EXP_D}/tleap/{vname}.prmtop"
inpcrd = f"{EXP_D}/tleap/{vname}.inpcrd"
chk = f"{md_dir}/{vname}_checkpoint.chk"

amber = pmd.load_file(prmtop, inpcrd)
system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                             constraints=app.HBonds, rigidWater=True)
system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))

ca_force = mm.CustomExternalForce("5.0 * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
ca_force.addPerParticleParameter("x0"); ca_force.addPerParticleParameter("y0"); ca_force.addPerParticleParameter("z0")
for atom in amber.atoms:
    if atom.residue.name not in ('WAT','HOH','SOL','Na+','Cl-','LNK') and atom.name == 'CA':
        xyz = amber.positions[atom.idx]; ca_force.addParticle(atom.idx, [xyz[0], xyz[1], xyz[2]])
system.addForce(ca_force)

integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
integrator.setRandomNumberSeed(rep * 42)
platform = mm.Platform.getPlatformByName('CUDA')
sim = app.Simulation(amber.topology, system, integrator, platform,
                      {'CudaDeviceIndex': args.gpu, 'CudaPrecision': 'mixed'})
sim.context.setPositions(amber.positions)

with open(chk, 'rb') as f:
    sim.context.loadCheckpoint(f.read())
state = sim.context.getState(getEnergy=True)
pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"Restored checkpoint. PE: {pe:.0f} kJ/mol", flush=True)

sim.reporters.append(app.DCDReporter(f"{md_dir}/{vname}_traj.dcd", 50000, append=True))
sim.reporters.append(app.StateDataReporter(f"{md_dir}/{vname}_log.txt", 10000,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, append=True))
sim.reporters.append(app.CheckpointReporter(f"{md_dir}/{vname}_checkpoint.chk", 500000))

production_steps = args.nsteps
sim_start = time.time(); steps_done = 0
while steps_done < production_steps:
    chunk = min(500000, production_steps - steps_done)
    sim.step(chunk); steps_done += chunk
    elapsed = time.time() - sim_start
    ns_done = steps_done * 2e-6; ns_day = ns_done / (elapsed / 86400)
    remaining = (production_steps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
    print(f"[{time.strftime('%H:%M:%S')}] +{ns_done:.0f}ns ({ns_day:.0f} ns/d, ~{remaining:.0f}h)", flush=True)
print("Done!")
