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

print("Minimizing...")
simulation.minimizeEnergy(maxIterations=100)
print("Minimization done")

print("Running 1 ps test...")
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(500)
print("Test done")
