#!/usr/bin/env python3
"""
Production MD for membrane-embedded GLP-1R + semaglutide (exp-F) - v2 with append support.
Restarts from production checkpoint, preserves existing DCD/log.

Input:  membrane_build/system_final.prmtop
        membrane_rep1/prod_checkpoint.chk (or membrane_equil/equil_safe.chk)
Output: membrane_rep1/prod_v2.{dcd,log,chk}
GPU:    CUDA device 1
"""
import sys, os, time, shutil
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
INPCRD = f"{EXP_F}/membrane_build/system_final.inpcrd"
OUT_DIR = f"{EXP_F}/md/membrane_rep1"
os.makedirs(OUT_DIR, exist_ok=True)

CHK_EQUIL = f"{EXP_F}/md/membrane_equil/equil_safe.chk"
CHK_PROD = f"{OUT_DIR}/prod_checkpoint.chk"

if os.path.exists(CHK_PROD) and os.path.exists(CHK_EQUIL):
    if os.path.getmtime(CHK_PROD) > os.path.getmtime(CHK_EQUIL):
        CHK_IN = CHK_PROD
    else:
        CHK_IN = CHK_EQUIL
elif os.path.exists(CHK_PROD):
    CHK_IN = CHK_PROD
else:
    CHK_IN = CHK_EQUIL

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
print("Membrane System Production MD v2")
print(f"Prmtop: {PRMTOP}")
print(f"Steps:  {N_STEPS} ({N_STEPS * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond):.1f} ns)")
print(f"Output: {OUT_DIR}")
print(f"GPU:    CUDA device 1")
print(f"Checkpoint: {CHK_IN}")
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

# Check current step
current_step = simulation.context.getState().getStepCount()
print(f"Current step from checkpoint: {current_step}")
remaining_steps = max(0, N_STEPS - current_step)
print(f"Remaining steps: {remaining_steps}")

if remaining_steps == 0:
    print("Production already complete!")
    sys.exit(0)

# Backup existing files if they exist and are large
for fname in ["prod.dcd", "prod.log"]:
    fpath = os.path.join(OUT_DIR, fname)
    if os.path.exists(fpath) and os.path.getsize(fpath) > 100_000_000:  # > 100 MB
        backup = fpath + ".backup_" + time.strftime("%Y%m%d_%H%M%S")
        shutil.copy2(fpath, backup)
        print(f"Backed up {fname} to {backup}")

# Reporters - write to NEW files to avoid overwriting
dcd_path = f"{OUT_DIR}/prod_v2.dcd"
log_path = f"{OUT_DIR}/prod_v2.log"
chk_path = f"{OUT_DIR}/prod_v2_checkpoint.chk"

simulation.reporters.append(app.DCDReporter(dcd_path, DCD_INTERVAL))
simulation.reporters.append(app.StateDataReporter(
    log_path, REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

# Checkpoint reporter
def write_checkpoint(step):
    with open(chk_path, "wb") as f:
        f.write(simulation.context.createCheckpoint())
    print(f"Checkpoint saved at step {step}")

print(f"\nStarting production: {remaining_steps} steps remaining")
print(f"Output DCD: {dcd_path}")
print(f"Output log: {log_path}")

t0 = time.time()

# Run in segments to allow checkpointing
segment = CHECKPOINT_INTERVAL
completed = current_step
for start in range(current_step, N_STEPS, segment):
    end = min(start + segment, N_STEPS)
    steps = end - start
    simulation.step(steps)
    completed = end
    write_checkpoint(end)
    elapsed = time.time() - t0
    speed_ns_day = ((completed - current_step) * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond)) / (elapsed / 86400)
    print(f"Completed {completed}/{N_STEPS} steps ({completed/N_STEPS*100:.1f}%) | Speed: {speed_ns_day:.1f} ns/day | Elapsed: {elapsed/3600:.1f} h")

print("\nProduction run complete!")
