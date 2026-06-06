#!/usr/bin/env python3
"""Minimization + short NVT MD for exp-F Phase 0 validation."""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
PRMTOP = f"{REPO}/exps/exp-F/structures/complex.prmtop"
INPCRD = f"{REPO}/exps/exp-F/structures/complex.inpcrd"
OUT = f"{REPO}/exps/exp-F/md"
os.makedirs(OUT, exist_ok=True)

# Load
print("Loading prmtop/inpcrd...")
amber = app.AmberPrmtopFile(PRMTOP)
inpcrd = app.AmberInpcrdFile(INPCRD)

# Create system
print("Creating system...")
system = amber.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True,
)

# Add weak restraints to receptor CA to prevent major drift during minimization
receptor_ca = []
for atom in amber.topology.atoms():
    if atom.name == 'CA' and atom.residue.chain.id in ('R', 'A', 'B', 'G', 'N'):
        receptor_ca.append(atom.index)

if receptor_ca:
    force = mm.CustomExternalForce("0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    force.addGlobalParameter("k", 1.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    for idx in receptor_ca:
        pos = inpcrd.positions[idx]
        force.addParticle(idx, [pos.x, pos.y, pos.z])
    system.addForce(force)
    print(f"Added weak restraints to {len(receptor_ca)} receptor CA atoms")

# Platform
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

# Minimization
integrator = mm.LangevinIntegrator(100*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
simulation = app.Simulation(amber.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

state = simulation.context.getState(getEnergy=True)
print(f"Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

print("Minimizing...")
simulation.minimizeEnergy(maxIterations=2000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"After minimization PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

# Save minimized structure
positions = state.getPositions()
with open(f"{OUT}/complex_minimized.pdb", "w") as f:
    app.PDBFile.writeFile(amber.topology, positions, f)
print(f"Saved minimized structure: {OUT}/complex_minimized.pdb")

# Short NVT heating + equilibration (10 ps at 100K, then 100 ps at 310K)
print("\nHeating 100K -> 310K (NVT, 110 ps total)...")
integrator.setTemperature(100*unit.kelvin)
simulation.step(5000)  # 10 ps at 100K

for temp in [150, 200, 250, 310]:
    integrator.setTemperature(temp*unit.kelvin)
    simulation.step(5000)  # 10 ps each step
    state = simulation.context.getState(getEnergy=True, getTemperature=True)
    print(f"  T={state.getTemperature().value_in_unit(unit.kelvin):.1f}K, PE={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

# Production: 1 ns NVT
print("\nProduction NVT: 1 ns")
report_interval = 5000  # every 10 ps
dcd_interval = 25000    # every 50 ps
checkpoint_interval = 500000  # every 1 ns

simulation.reporters.append(app.DCDReporter(f"{OUT}/complex_traj.dcd", dcd_interval))
simulation.reporters.append(app.StateDataReporter(
    f"{OUT}/complex_log.txt", report_interval,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    temperature=True, volume=True, density=True, speed=True
))
simulation.reporters.append(app.CheckpointReporter(f"{OUT}/complex_checkpoint.chk", checkpoint_interval))

nsteps = 500000  # 1 ns at 2 fs
sim_start = time.time()
simulation.step(nsteps)

elapsed = time.time() - sim_start
ns_done = nsteps * 2e-6
ns_day = ns_done / (elapsed / 86400)
print(f"\nDone! {ns_done:.0f} ns in {elapsed/60:.1f} min ({ns_day:.0f} ns/day)")
