#!/usr/bin/env python3
"""Safe membrane equilibration with exception handling."""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
INPCRD = f"{EXP_F}/membrane_build/system_final.inpcrd"
OUT_DIR = f"{EXP_F}/md/membrane_equil"
os.makedirs(OUT_DIR, exist_ok=True)

platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "1"}

print("Loading...")
prmtop = app.AmberPrmtopFile(PRMTOP)
inpcrd = app.AmberInpcrdFile(INPCRD)

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)

# CA restraints
receptor_ca = [a.index for a in prmtop.topology.atoms() if a.name == 'CA' and a.residue.index < 1113]
if receptor_ca:
    f = mm.CustomExternalForce("0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2")
    f.addPerParticleParameter("x0"); f.addPerParticleParameter("y0"); f.addPerParticleParameter("z0")
    f.addGlobalParameter("k", 10.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    for idx in receptor_ca:
        p = inpcrd.positions[idx]
        f.addParticle(idx, [p.x, p.y, p.z])
    system.addForce(f)
    print(f"Restraints on {len(receptor_ca)} CA")

integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platform_props)
simulation.context.setPositions(inpcrd.positions)

# Load minimized checkpoint if exists
min_pdb = f"{OUT_DIR}/minimized.pdb"
if os.path.exists(min_pdb):
    print(f"Loading minimized positions from {min_pdb}")
    min_struct = app.PDBFile(min_pdb)
    simulation.context.setPositions(min_struct.positions)

# Add barostat
system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 310 * unit.kelvin))
simulation.context.reinitialize(preserveState=True)

simulation.reporters.append(app.StateDataReporter(
    f"{OUT_DIR}/equil_safe.log", 5000,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

print("Running NPT equilibration (100 ps)...")
try:
    simulation.step(50000)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"Final PE: {state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):,.2f}")
    with open(f"{OUT_DIR}/equilibrated_safe.pdb", "w") as f:
        app.PDBFile.writeFile(prmtop.topology, state.getPositions(), f)
    with open(f"{OUT_DIR}/equil_safe.chk", "wb") as f:
        f.write(simulation.context.createCheckpoint())
    print("Success!")
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
