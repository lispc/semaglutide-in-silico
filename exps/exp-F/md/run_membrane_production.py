#!/usr/bin/env python3
"""
Production MD for membrane-embedded GLP-1R + semaglutide (exp-F).
Restarts from membrane_equil checkpoint.

Input:  membrane_build/system_final.prmtop
        membrane_equil/equil_safe.chk
Output: membrane_rep1/prod.{dcd,log,chk}
GPU:    CUDA device 1
"""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
INPCRD = f"{EXP_F}/membrane_build/system_final.inpcrd"
CHK_IN = f"{EXP_F}/md/membrane_equil/equil_safe.chk"
OUT_DIR = f"{EXP_F}/md/membrane_rep1"
os.makedirs(OUT_DIR, exist_ok=True)

# Simulation parameters
TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond
N_STEPS = 50_000_000  # 100 ns
REPORT_INTERVAL = 5000  # 10 ps
CHECKPOINT_INTERVAL = 500_000  # 1 ns
DCD_INTERVAL = 5000  # 10 ps

# Platform (GPU 1, leave GPU 0 for solvent production)
platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "1"}

print("=" * 60)
print("Membrane System Production MD")
print(f"Prmtop: {PRMTOP}")
print(f"Steps:  {N_STEPS} ({N_STEPS * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond):.1f} ns)")
print(f"Output: {OUT_DIR}")
print(f"GPU:    CUDA device 1")
print("=" * 60)

print("\nLoading system...")
prmtop = app.AmberPrmtopFile(PRMTOP)
inpcrd = app.AmberInpcrdFile(INPCRD)

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005,
)

# Optional: weak restraints on receptor CA during production
# Set to False to run unrestrained production
USE_CA_RESTRAINTS = False

if USE_CA_RESTRAINTS:
    receptor_ca = []
    for atom in prmtop.topology.atoms():
        if atom.name == 'CA' and atom.residue.chain.id in ('R', 'A', 'B', 'G', 'N'):
            receptor_ca.append(atom.index)
    if receptor_ca:
        restraint_force = mm.CustomExternalForce(
            "0.5 * k * periodicdistance(x, y, z, x0, y0, z0)^2"
        )
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        restraint_force.addGlobalParameter("k", 1.0 * unit.kilojoule_per_mole / unit.nanometer**2)
        for idx in receptor_ca:
            pos = inpcrd.positions[idx]
            restraint_force.addParticle(idx, [pos.x, pos.y, pos.z])
        system.addForce(restraint_force)
        print(f"Added weak CA restraints ({len(receptor_ca)} atoms)")

# NPT
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platform_props)
simulation.context.setPositions(inpcrd.positions)

# Load checkpoint
if os.path.exists(CHK_IN):
    with open(CHK_IN, "rb") as f:
        simulation.context.loadCheckpoint(f.read())
    print(f"Loaded checkpoint from {CHK_IN}")
else:
    print(f"WARNING: No checkpoint found at {CHK_IN}, starting from inpcrd")
    sys.exit(1)

# Reporters
dcd_path = f"{OUT_DIR}/prod.dcd"
log_path = f"{OUT_DIR}/prod.log"
chk_path = f"{OUT_DIR}/prod_checkpoint.chk"

simulation.reporters.append(app.DCDReporter(dcd_path, DCD_INTERVAL))
simulation.reporters.append(app.StateDataReporter(
    log_path, REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

def write_checkpoint(step):
    with open(chk_path, "wb") as f:
        f.write(simulation.context.createCheckpoint())
    print(f"Checkpoint saved at step {step}")

print(f"\nStarting production: {N_STEPS} steps")
print(f"Output DCD: {dcd_path}")
print(f"Output log: {log_path}")

t0 = time.time()

# Run in segments to allow checkpointing
segment = CHECKPOINT_INTERVAL
for start in range(0, N_STEPS, segment):
    end = min(start + segment, N_STEPS)
    steps = end - start
    simulation.step(steps)
    write_checkpoint(end)
    elapsed = time.time() - t0
    speed_ns_day = (end * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond)) / (elapsed / 86400)
    print(f"Completed {end}/{N_STEPS} steps ({end/N_STEPS*100:.1f}%) | Speed: {speed_ns_day:.1f} ns/day | Elapsed: {elapsed/3600:.1f} h")

print("\nProduction run complete!")
