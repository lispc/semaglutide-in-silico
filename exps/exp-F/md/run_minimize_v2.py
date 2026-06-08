#!/usr/bin/env python3
"""Minimization + NVT heating/equilibration for ECD v2 model."""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
PRMTOP = f"{REPO}/exps/exp-F/minimal_model/system_ecd_v2.prmtop"
INPCRD = f"{REPO}/exps/exp-F/minimal_model/system_ecd_v2.inpcrd"
OUT = f"{REPO}/exps/exp-F/md/ecd_v2"
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

# Add weak restraints to receptor and HSA CA to prevent major drift during minimization
receptor_ca = []
for atom in amber.topology.atoms():
    if atom.name == 'CA' and atom.residue.chain.id in ('R', 'A'):
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
    print(f"Added weak restraints to {len(receptor_ca)} CA atoms")

# Platform
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

# Minimization
integrator = mm.LangevinIntegrator(100*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
simulation = app.Simulation(amber.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)

state = simulation.context.getState(getEnergy=True)
print(f"Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

print("Minimizing (max 5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"After minimization PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

# Save minimized structure
positions = state.getPositions()
with open(f"{OUT}/minimized.pdb", "w") as f:
    app.PDBFile.writeFile(amber.topology, positions, f)
print(f"Saved minimized structure: {OUT}/minimized.pdb")

# Save checkpoint
with open(f"{OUT}/minimized.chk", "wb") as f:
    f.write(simulation.context.createCheckpoint())
print(f"Saved checkpoint: {OUT}/minimized.chk")

# Short NVT heating (100K -> 310K, 100 ps total)
print("\nHeating 100K -> 310K (NVT, 100 ps)...")
integrator.setTemperature(100*unit.kelvin)
for step, temp in enumerate([100, 150, 200, 250, 310]):
    integrator.setTemperature(temp*unit.kelvin)
    simulation.step(10000)  # 20 ps each step
    state = simulation.context.getState(getEnergy=True)
    # Approximate temperature from KE
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    T = (2 * state.getKineticEnergy()) / (dof * unit.MOLAR_GAS_CONSTANT_R)
    print(f"  Step {step+1}/5: T={T.value_in_unit(unit.kelvin):.1f}K, PE={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.2f} kJ/mol")

# Short NVT equilibration (1 ns)
print("\nNVT equilibration: 1 ns")
simulation.reporters.append(app.DCDReporter(f"{OUT}/equil.dcd", 5000))  # 10 ps/frame
simulation.reporters.append(app.StateDataReporter(
    f"{OUT}/equil.log", 5000,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    temperature=True, volume=True, density=True, speed=True
))
simulation.step(500000)  # 1 ns

# Save equilibration checkpoint
with open(f"{OUT}/equil.chk", "wb") as f:
    f.write(simulation.context.createCheckpoint())
print(f"Saved equil checkpoint: {OUT}/equil.chk")

print("\nMinimization and equilibration complete!")
