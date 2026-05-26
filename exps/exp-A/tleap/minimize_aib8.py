#!/usr/bin/env python3
"""Minimize Aib8 system (ParmEd modified topology) and prepare for production."""
import parmed as pmd, openmm as mm, openmm.app as app, openmm.unit as unit, numpy as np, os

out_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/md/aib8"
amber = pmd.load_file(f"{out_dir}/aib8_modified.prmtop", f"{out_dir}/aib8_modified.inpcrd")

# Create OpenMM system via ParmEd (has all params including AIB)
system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                             constraints=app.HBonds, rigidWater=True)
with open(f"{out_dir}/aib8_system.xml", "w") as f:
    f.write(mm.XmlSerializer.serialize(system))

# Save PDB
amber.save(f"{out_dir}/aib8_amber.pdb", overwrite=True)
pdb = app.PDBFile(f"{out_dir}/aib8_amber.pdb")

# Minimize
integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
sim = app.Simulation(pdb.topology, system, integrator)
sim.context.setPositions(amber.positions)

pe0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
print(f"PE before min: {pe0}")
sim.minimizeEnergy(maxIterations=5000)
pe1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
print(f"PE after min: {pe1}")

state = sim.context.getState(getPositions=True)
app.PDBFile.writeFile(pdb.topology, state.getPositions(), open(f"{out_dir}/aib8_minimized.pdb", "w"))
print(f"Minimized saved: {out_dir}/aib8_minimized.pdb")
print("Done!")
