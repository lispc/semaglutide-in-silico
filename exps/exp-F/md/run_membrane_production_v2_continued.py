#!/usr/bin/env python3
"""
Production MD continuation for membrane-embedded GLP-1R + semaglutide (exp-F).
Uses OLD system (system_old_v3.prmtop) to continue from 110 ns to 200 ns.
The ONLY difference vs current system is 1 DUM residue (21 virtual atoms).

Input:  temp_old_system/system_old_v3.prmtop
        temp_old_system/system_old_v3.inpcrd (or old_defective/prod_v2_checkpoint.chk)
Output: membrane_rep1/prod_v2_cont.{dcd,log,chk}
GPU:    CUDA device 1
"""
import sys, os, time, shutil
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/temp_old_system/system_old_v3.prmtop"
INPCRD = f"{EXP_F}/temp_old_system/system_old.inpcrd"  # inpcrd for topology init; actual coords come from checkpoint
OUT_DIR = f"{EXP_F}/md/membrane_rep1"
os.makedirs(OUT_DIR, exist_ok=True)

# Checkpoint sources (priority: current dir > old_defective)
CHK_NEW = f"{OUT_DIR}/prod_v2_cont_checkpoint.chk"
CHK_OLD = f"{OUT_DIR}/old_defective/prod_v2_checkpoint.chk"

if os.path.exists(CHK_NEW) and os.path.getsize(CHK_NEW) > 1000000:
    CHK_IN = CHK_NEW
    print(f"Using current continuation checkpoint: {CHK_IN}")
elif os.path.exists(CHK_OLD) and os.path.getsize(CHK_OLD) > 1000000:
    CHK_IN = CHK_OLD
    print(f"Using old 110ns checkpoint: {CHK_IN}")
else:
    print("ERROR: No valid checkpoint found!")
    sys.exit(1)

# Simulation parameters
TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond
N_STEPS = 200_000_000  # 200 ns target
REPORT_INTERVAL = 5000  # 10 ps
CHECKPOINT_INTERVAL = 500_000  # 1 ns
DCD_INTERVAL = 5000  # 10 ps

platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "1"}

print("=" * 60)
print("Membrane Production MD - CONTINUATION (110 ns -> 200 ns)")
print(f"Prmtop: {PRMTOP}")
print(f"Steps target: {N_STEPS} (200 ns)")
print(f"Checkpoint: {CHK_IN}")
print("NOTE: Using old system with DUM residue (21 virtual atoms)")
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
system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE))

integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
simulation = app.Simulation(prmtop.topology, system, integrator, platform, platform_props)
# Note: inpcrd has fewer atoms (no DUM), so skip setPositions.
# Checkpoint contains full coordinates including DUM atoms.
# with open(CHK_IN, "rb") as f:
#     simulation.context.loadCheckpoint(f.read())

# Load checkpoint
with open(CHK_IN, "rb") as f:
    simulation.context.loadCheckpoint(f.read())
print(f"Loaded checkpoint from {CHK_IN}")

current_step = simulation.context.getState().getStepCount()
print(f"Current step: {current_step:,} -> {current_step * 0.002 / 1000:.2f} ns")
remaining_steps = max(0, N_STEPS - current_step)
print(f"Remaining steps: {remaining_steps:,} -> {remaining_steps * 0.002 / 1000:.2f} ns")

if remaining_steps == 0:
    print("Production already complete!")
    sys.exit(0)

# Output to NEW files (prod_v2_cont) to avoid overwriting anything
dcd_path = f"{OUT_DIR}/prod_v2_cont.dcd"
log_path = f"{OUT_DIR}/prod_v2_cont.log"
chk_path = f"{OUT_DIR}/prod_v2_cont_checkpoint.chk"

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

print(f"\nStarting production: {remaining_steps} steps remaining")
print(f"Output DCD: {dcd_path}")
print(f"Output log: {log_path}")

t0 = time.time()
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
