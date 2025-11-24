import MDAnalysis as mda
import numpy as np
from itertools import combinations

# --- INPUT FILES ---
TOPOLOGY = "topology.pdb"
TRAJECTORY = "traj.xtc"

u = mda.Universe(TOPOLOGY, TRAJECTORY)

# --- HELIX DEFINITIONS (Cα atoms only) ---
helix_resid_ranges = [
    (22, 33),     # Helix 1
    (286, 297),   # Helix 2
    (350, 355),   # Helix 3
]

helices = []
for start, end in helix_resid_ranges:
    sel = u.select_atoms(f"resid {start}-{end} and name CA")
    if len(sel) == 0:
        raise ValueError(f"No atoms found for resid {start}-{end}")
    helices.append(sel)

# All helix pairs
pairs = list(combinations(range(len(helices)), 2))

# Arrays to store results
times = []
mean_pair_dists = []

# --- LOOP THROUGH TRAJECTORY ---
for ts in u.trajectory:
    times.append(ts.time)

    # Compute COM of each helix
    coms = np.array([h.center_of_mass() for h in helices])

    # Compute pairwise distances between COMs
    dists = []
    for i, j in pairs:
        d = np.linalg.norm(coms[i] - coms[j])
        dists.append(d)
    dists = np.array(dists)

    # Mean pairwise distance = global openness
    mean_pair_dists.append(dists.mean())

mean_pair_dists = np.array(mean_pair_dists)
times = np.array(times)

# --- FIND THE MOST OPEN FRAME ---
frame_open = int(np.argmax(mean_pair_dists))

print("=== MOST OPEN CRYPTIC-POCKET FRAME (max mean COM distance) ===")
print(f"Frame index: {frame_open}")
print(f"Time: {times[frame_open]:.3f}")
print(f"Mean helix COM distance: {mean_pair_dists[frame_open]:.3f} Å")

# --- OPTIONAL: WRITE OUT THE FRAME FOR VISUALIZATION ---
u.trajectory[frame_open]
u.atoms.write("most_open_cryptic_pocket.pdb")
print("\nSaved structure: most_open_cryptic_pocket.pdb")
