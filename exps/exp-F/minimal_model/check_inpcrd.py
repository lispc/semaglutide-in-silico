from openmm.app import AmberInpcrdFile
import numpy as np

inpcrd = AmberInpcrdFile('system_ecd.inpcrd')
pos = inpcrd.positions
print(f"Number of positions: {len(pos)}")
print(f"First position: {pos[0]}")
print(f"Last position: {pos[-1]}")

# Convert to Angstrom manually
vals = np.array([[p.x, p.y, p.z] for p in pos])
# OpenMM uses nanometers internally, but AmberInpcrdFile may return angstrom
# Let's check by looking at the raw values
print(f"X range: {vals[:,0].min():.4f} to {vals[:,0].max():.4f}")
print(f"Y range: {vals[:,1].min():.4f} to {vals[:,1].max():.4f}")
print(f"Z range: {vals[:,2].min():.4f} to {vals[:,2].max():.4f}")

# If max is ~167, unit is angstrom (box Z is 167 A)
# If max is ~16.7, unit is nanometer
print(f"Box vectors: {inpcrd.boxVectors}")
