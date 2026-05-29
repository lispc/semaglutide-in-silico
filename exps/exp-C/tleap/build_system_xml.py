#!/usr/bin/env python3
"""Build and cache OpenMM system XML to avoid repeated createSystem calls."""
import parmed as pmd, openmm.app as app, openmm.unit as unit, openmm as mm, sys, os

sys_name = sys.argv[1] if len(sys.argv) > 1 else "linker_c18"
TL_DIR = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap"
prmtop = f"{TL_DIR}/{sys_name}.prmtop"
inpcrd = f"{TL_DIR}/{sys_name}.inpcrd"
system_xml = f"{TL_DIR}/{sys_name}_system.xml"

if os.path.exists(system_xml):
    print(f"{system_xml} already exists, skipping")
    sys.exit(0)

print(f"Loading {prmtop}...")
amber = pmd.load_file(prmtop, inpcrd)
print(f"Building system for {len(amber.atoms)} atoms (this may take 1-2 min)...")
system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                             constraints=app.HBonds, rigidWater=True)
print(f"Saving to {system_xml}...")
with open(system_xml, 'w') as f:
    f.write(mm.XmlSerializer.serialize(system))
print("Done!")
