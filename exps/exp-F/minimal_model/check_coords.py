from openmm.app import AmberInpcrdFile
from openmm.unit import angstrom, nanometer
import numpy as np

inpcrd = AmberInpcrdFile('system_ecd.inpcrd')
pos = np.array([[p.x, p.y, p.z] for p in inpcrd.positions])  # in nm
pos_ang = pos * 10
print(f"Shape: {pos_ang.shape}")
print(f"X range: {pos_ang[:,0].min():.2f} to {pos_ang[:,0].max():.2f}")
print(f"Y range: {pos_ang[:,1].min():.2f} to {pos_ang[:,1].max():.2f}")
print(f"Z range: {pos_ang[:,2].min():.2f} to {pos_ang[:,2].max():.2f}")
print(f"Any NaN: {np.any(np.isnan(pos_ang))}")

# Check pairwise distances for a sample
import random
n = len(pos_ang)
sample_size = min(2000, n)
idx = random.sample(range(n), sample_size)
sample = pos_ang[idx]

# Compute min distance for each point
from scipy.spatial import cKDTree
tree = cKDTree(sample)
distances, _ = tree.query(sample, k=2)  # k=2 because nearest neighbor includes self
min_distances = distances[:, 1]
print(f"Sample size: {sample_size}")
print(f"Min distance in sample: {min_distances.min():.4f} A")
print(f"Atoms with dist < 0.5 A: {(min_distances < 0.5).sum()}")
print(f"Atoms with dist < 1.0 A: {(min_distances < 1.0).sum()}")

# Check box vectors
if inpcrd.boxVectors is not None:
    bv = inpcrd.boxVectors
    print(f"Box vectors: {[v.value_in_unit(angstrom) for v in bv]}")
