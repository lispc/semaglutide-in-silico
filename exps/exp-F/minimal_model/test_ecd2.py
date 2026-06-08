from openmm import *
from openmm.app import *
from openmm.unit import *
import sys

prmtop = AmberPrmtopFile('system_ecd.prmtop')
inpcrd = AmberInpcrdFile('system_ecd.inpcrd')

system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True,
)

system.addForce(MonteCarloBarostat(1.0*bar, 300*kelvin))
integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 0.002*picosecond)
platform = Platform.getPlatformByName('CUDA')
simulation = Simulation(prmtop.topology, system, integrator, platform, {'CudaDeviceIndex': '2', 'CudaPrecision': 'mixed'})

simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Check initial energy
state = simulation.context.getState(getEnergy=True)
print(f"Initial energy: {state.getPotentialEnergy()}")

print("Minimizing (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000, tolerance=10*kilojoule_per_mole)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"Minimized energy: {state.getPotentialEnergy()}")

# Check for NaN
pos = state.getPositions(asNumpy=True)
import numpy as np
if np.any(np.isnan(pos.value_in_unit(angstrom))):
    print("WARNING: NaN coordinates found after minimization!")
    nan_mask = np.isnan(pos.value_in_unit(angstrom))
    nan_atoms = np.where(np.any(nan_mask, axis=1))[0]
    print(f"NaN atoms: {nan_atoms[:10]} (total {len(nan_atoms)})")
else:
    print("No NaN coordinates after minimization.")

# Save minimized structure
with open('minimized_test.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

print("Running 10 ps test...")
simulation.context.setVelocitiesToTemperature(300*kelvin)
for i in range(10):
    simulation.step(500)  # 1 ps
    state = simulation.context.getState(getEnergy=True)
    print(f"Step {(i+1)*500}: PE={state.getPotentialEnergy()}")
print("Test done")
