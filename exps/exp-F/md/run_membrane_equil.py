#!/usr/bin/env python3
"""
Membrane system equilibration for exp-F.
Minimization -> NVT heating -> NPT heating -> NPT equilibration.
Restraints on receptor CA during heating/equil.

Input:  membrane_build/system_final.{prmtop,inpcrd}
Output: membrane_equil/{minimized.pdb,equilibrated.pdb,equil.chk,equil.log}
GPU:    CUDA device 1
"""
import sys, os, time
start_time = time.time()
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
INPCRD = f"{EXP_F}/membrane_build/system_final.inpcrd"
OUT_DIR = f"{EXP_F}/md/membrane_equil"
os.makedirs(OUT_DIR, exist_ok=True)

# Simulation parameters
TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond

# Platform (use GPU 1, leave GPU 0 for solvent production)
platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "1"}

print("=" * 60)
print("Membrane System Equilibration")
print(f"Prmtop: {PRMTOP}")
print(f"Inpcrd: {INPCRD}")
print(f"Output: {OUT_DIR}")
print(f"GPU: CUDA device 1")
print("=" * 60)

# Load topology and coordinates
print("\n[1/5] Loading system...")
prmtop = app.AmberPrmtopFile(PRMTOP)
inpcrd = app.AmberInpcrdFile(INPCRD)

# Create system
print("[2/5] Creating OpenMM system...")
system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)

# Identify receptor CA atoms for restraints
# tleap merges everything into one chain; use residue index to distinguish
# Receptor: residues 1-1113; Peptide (semaglutide): residues 1114-1140
receptor_ca = []
for atom in prmtop.topology.atoms():
    if atom.name == 'CA' and atom.residue.index < 1113:
        receptor_ca.append(atom.index)

print(f"  Receptor CA atoms for restraint: {len(receptor_ca)}")

# Add weak harmonic restraints to receptor CA
if receptor_ca:
    restraint_force = mm.CustomExternalForce(
        "0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2"
    )
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    restraint_force.addGlobalParameter("k", 10.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    for idx in receptor_ca:
        pos = inpcrd.positions[idx]
        restraint_force.addParticle(idx, [pos.x, pos.y, pos.z])
    system.addForce(restraint_force)
    print("  Added 10 kJ/mol/nm^2 restraints to receptor CA")

# Minimization (no barostat)
print("\n[3/5] Energy minimization...")
integrator = mm.LangevinMiddleIntegrator(100 * unit.kelvin, FRICTION, TIMESTEP)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platform_props)
simulation.context.setPositions(inpcrd.positions)

state = simulation.context.getState(getEnergy=True)
print(f"  Initial PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f} kJ/mol")

simulation.minimizeEnergy(maxIterations=1000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"  After minimization PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f} kJ/mol")

# Save minimized structure
positions = state.getPositions()
with open(f"{OUT_DIR}/minimized.pdb", "w") as f:
    app.PDBFile.writeFile(prmtop.topology, positions, f)
print(f"  Saved: {OUT_DIR}/minimized.pdb")

# NVT heating: 0 -> 100 K, 50 ps
print("\n[4/5] NVT heating (0 -> 100 K, 50 ps)...")
integrator.setTemperature(100 * unit.kelvin)
report_interval = 5000  # 10 ps

simulation.reporters.append(app.StateDataReporter(
    f"{OUT_DIR}/heat_nvt.log", report_interval,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

simulation.step(25000)  # 50 ps
state = simulation.context.getState(getEnergy=True)
print(f"  PE={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f} kJ/mol")

# NPT heating + equilibration
print("\n[5/5] NPT heating (100 -> 310 K) + equilibration...")
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE))
simulation.context.reinitialize(preserveState=True)

# Remove old reporter and add new one
simulation.reporters.clear()
simulation.reporters.append(app.StateDataReporter(
    f"{OUT_DIR}/equil.log", report_interval,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

# Heating stages: 100 -> 150 -> 200 -> 250 -> 310 K
heat_temps = [150, 200, 250, 310]
steps_per_stage = 10000  # 20 ps each -> 80 ps heating

for temp in heat_temps:
    integrator.setTemperature(temp * unit.kelvin)
    simulation.context.setParameter(mm.MonteCarloBarostat.Temperature(), temp * unit.kelvin)
    simulation.step(steps_per_stage)
    state = simulation.context.getState(getEnergy=True)
    print(f"  PE={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f} kJ/mol, "
          f"V={state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3):.1f} nm^3")

# NPT equilibration at 310 K
print("\n  NPT equilibration at 310 K (100 ps)...")
simulation.step(50000)  # 100 ps

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"  Final: PE={state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f} kJ/mol, "
      f"V={state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3):.1f} nm^3")

# Save equilibrated structure and checkpoint
positions = state.getPositions()
with open(f"{OUT_DIR}/equilibrated.pdb", "w") as f:
    app.PDBFile.writeFile(prmtop.topology, positions, f)

chk_path = f"{OUT_DIR}/equil.chk"
with open(chk_path, "wb") as f:
    f.write(simulation.context.createCheckpoint())
print(f"  Saved: {OUT_DIR}/equilibrated.pdb")
print(f"  Saved checkpoint: {chk_path}")

elapsed = time.time() - start_time
print("\n" + "=" * 60)
print("Equilibration complete!")
print(f"Total time: {elapsed/60:.1f} min")
print("Next: run_membrane_production.py")
print("=" * 60)
