#!/usr/bin/env python3
"""Pre-flight OpenMM minimization check for exp-D rebuilt systems (2026-07-17).
Loads tleap/{v}.prmtop/.inpcrd, minimizes 500 steps on CPU, reports PE
before/after. Expect ~-5.7e5..-6.0e5 kJ/mol after min (exp-log reference),
anything >= 0 or NaN indicates a broken start.
"""
import sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"

v = sys.argv[1]
amber = pmd.load_file(f"{TLEAP}/{v}.prmtop", f"{TLEAP}/{v}.inpcrd")
system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometers,
                            constraints=app.HBonds, rigidWater=True)
integ = mm.LangevinIntegrator(310 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds)
sim = app.Simulation(amber.topology, system, integ, mm.Platform.getPlatformByName("CPU"))
sim.context.setPositions(amber.positions)
st = sim.context.getState(getEnergy=True)
pe0 = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
sim.minimizeEnergy(maxIterations=500)
st = sim.context.getState(getEnergy=True)
pe1 = st.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"{v}: PE init {pe0:.0f} -> PE min500 {pe1:.0f} kJ/mol", flush=True)
