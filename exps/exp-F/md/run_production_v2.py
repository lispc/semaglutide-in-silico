#!/usr/bin/env python3
"""
Production MD for ECD v2 model: 100 ns NPT at 310 K.
Restarts from equilibrated checkpoint.

Output:
  - md/ecd_v2/prod.dcd (trajectory, 10 ps/frame)
  - md/ecd_v2/prod.log (energy log, 10 ps)
  - md/ecd_v2/prod_checkpoint.chk (restart, every 1 ns)
"""
import sys, os, time
import openmm as mm
import openmm.app as app
import openmm.unit as unit

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/minimal_model/system_ecd_v2.prmtop"
INPCRD = f"{EXP_F}/minimal_model/system_ecd_v2.inpcrd"
CHK_IN = f"{EXP_F}/md/ecd_v2/equil.chk"
OUT_DIR = f"{EXP_F}/md/ecd_v2"

# Simulation parameters
TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond
N_STEPS = 50_000_000  # 100 ns
REPORT_INTERVAL = 5000  # 10 ps
CHECKPOINT_INTERVAL = 500_000  # 1 ns
DCD_INTERVAL = 5000  # 10 ps

# Platform
platform = mm.Platform.getPlatformByName("CUDA")
platform_props = {"CudaPrecision": "mixed", "CudaDeviceIndex": "0"}

print("Loading system...")
prmtop = app.AmberPrmtopFile(PRMTOP)
inpcrd = app.AmberInpcrdFile(INPCRD)

system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
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

# Reporters
dcd_path = f"{OUT_DIR}/prod.dcd"
log_path = f"{OUT_DIR}/prod.log"
chk_path = f"{OUT_DIR}/prod_checkpoint.chk"

simulation.reporters.append(app.DCDReporter(dcd_path, DCD_INTERVAL))
simulation.reporters.append(app.StateDataReporter(
    log_path,
    REPORT_INTERVAL,
    step=True,
    time=True,
    potentialEnergy=True,
    kineticEnergy=True,
    totalEnergy=True,
    temperature=True,
    volume=True,
    density=True,
    speed=True,
    separator=",",
))

# Checkpoint reporter
def write_checkpoint(step):
    with open(chk_path, "wb") as f:
        f.write(simulation.context.createCheckpoint())
    print(f"Checkpoint saved at step {step}")

print(f"Starting production run: {N_STEPS} steps ({N_STEPS * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond):.1f} ns)")
print(f"Output DCD: {dcd_path}")
print(f"Output log: {log_path}")

t0 = time.time()

# Run in segments to allow checkpointing
segment = CHECKPOINT_INTERVAL
for start in range(current_step, N_STEPS, segment):
    end = min(start + segment, N_STEPS)
    steps = end - start
    simulation.step(steps)
    write_checkpoint(end)
    elapsed = time.time() - t0
    speed_ns_day = (end * TIMESTEP.in_units_of(unit.nanosecond).value_in_unit(unit.nanosecond)) / (elapsed / 86400)
    print(f"Completed {end}/{N_STEPS} steps ({end/N_STEPS*100:.1f}%) | Speed: {speed_ns_day:.1f} ns/day | Elapsed: {elapsed/3600:.1f} h")

print("Production run complete!")
