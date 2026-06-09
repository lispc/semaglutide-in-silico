#!/usr/bin/env python3
"""
Extend membrane production MD from 100 ns -> 200 ns.
Restarts from prod_v2_checkpoint.chk, appends to existing DCD/log.

GPU: CUDA device 1
"""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
INPCRD = f"{EXP_F}/membrane_build/system_final.inpcrd"
OUT_DIR = f"{EXP_F}/md/membrane_rep1"
CHK_IN = f"{OUT_DIR}/prod_v2_checkpoint.chk"

# Extended to 200 ns
N_STEPS = 100_000_000  # 200 ns
REPORT_INTERVAL = 5000  # 10 ps
CHECKPOINT_INTERVAL = 500_000  # 1 ns
DCD_INTERVAL = 5000  # 10 ps

TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond

platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "1"}

print("=" * 60)
print("Membrane System Production MD EXTENSION 100 ns -> 200 ns")
print(f"Prmtop: {PRMTOP}")
print(f"Target: {N_STEPS} steps ({N_STEPS * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond):.1f} ns)")
print(f"GPU:    CUDA device 1")
print(f"Checkpoint: {CHK_IN}")
print("=" * 60)

if not os.path.exists(CHK_IN):
    print(f"ERROR: Checkpoint not found: {CHK_IN}")
    sys.exit(1)

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
simulation.context.setPositions(inpcrd.positions)

with open(CHK_IN, "rb") as f:
    simulation.context.loadCheckpoint(f.read())
print(f"Loaded checkpoint from {CHK_IN}")

current_step = simulation.context.getState().getStepCount()
print(f"Current step from checkpoint: {current_step}")
remaining_steps = max(0, N_STEPS - current_step)
print(f"Remaining steps: {remaining_steps}")

if remaining_steps == 0:
    print("Target already reached!")
    sys.exit(0)

# Append to existing DCD and log
dcd_path = f"{OUT_DIR}/prod_v2.dcd"
log_path = f"{OUT_DIR}/prod_v2.log"
chk_path = f"{OUT_DIR}/prod_v2_checkpoint.chk"

print(f"Appending to DCD: {dcd_path}")
print(f"Appending to log: {log_path}")

simulation.reporters.append(app.DCDReporter(dcd_path, DCD_INTERVAL, append=True))

# StateDataReporter with append mode (open file in append mode)
log_file = open(log_path, "a")
simulation.reporters.append(app.StateDataReporter(
    log_file, REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, kineticEnergy=True,
    totalEnergy=True, temperature=True, volume=True, density=True, speed=True,
    separator=",",
))

def write_checkpoint(step):
    with open(chk_path, "wb") as f:
        f.write(simulation.context.createCheckpoint())
    print(f"Checkpoint saved at step {step}")

print(f"\nStarting extension: {remaining_steps} steps remaining")
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

log_file.close()
print("\nExtension complete! Total: 200 ns")
