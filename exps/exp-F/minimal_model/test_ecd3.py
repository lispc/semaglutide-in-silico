from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
import numpy as np

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

state = simulation.context.getState(getEnergy=True)
print(f"Initial energy: {state.getPotentialEnergy()}")

print("Minimizing (5000 steps)...")
simulation.minimizeEnergy(maxIterations=5000)

state = simulation.context.getState(getEnergy=True, getPositions=True)
pe = state.getPotentialEnergy()
print(f"Minimized energy: {pe}")

pos = state.getPositions(asNumpy=True).value_in_unit(angstrom)
if np.any(np.isnan(pos)):
    print("WARNING: NaN found after minimization!")
else:
    print("No NaN after minimization.")

# Check if energy is reasonable
pe_kj = pe.value_in_unit(kilojoule_per_mole)
if pe_kj > 1e9 or np.isnan(pe_kj):
    print(f"Energy still too high or NaN: {pe_kj}. Not running MD.")
    sys.exit(1)

print("Running 10 ps test...")
simulation.context.setVelocitiesToTemperature(300*kelvin)
for i in range(10):
    simulation.step(500)
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy()
    pe_kj = pe.value_in_unit(kilojoule_per_mole)
    print(f"Step {(i+1)*500}: {pe}")
    if np.isnan(pe_kj):
        print("NaN energy during MD!")
        break

print("Test done")
