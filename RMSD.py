"""
Local Backbone RMSD Calculator using MDAnalysis.

Description:
    This script calculates the RMSD of protein backbone segments (3-residue sliding window)
    relative to a reference structure (topology file).
    It processes all matching trajectory files in the current directory,
    computes the RMSD frame-by-frame, and then bins the data into time windows.
    
    Output: Separate CSV files containing average RMSD per residue for each time window.

Requirements:
    - MDAnalysis
    - numpy

Usage:
    Place this script in the directory containing your .gro and .xtc files and run:
    python RMSD.py
"""

import os
import csv
import warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms

# Suppress specific MDAnalysis warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)

# --- CONFIGURATION ---
WINDOW_SIZE_PS = 50000.0       # Time window for averaging (50 ns = 50,000 ps)
TOPOLOGY_FILE = "step7_1.gro"  # Reference/Topology file
TRAJECTORY_PREFIX = "prod"     # Prefix for trajectory files (e.g., "prod")
TRAJECTORY_EXT = ".xtc"        # Extension for trajectory files
# ---------------------

def main():
    folder_path = os.getcwd()

    print(f"--- Local RMSD Analysis ---")
    print(f"Working directory: {folder_path}")
    print(f"Time window size:  {WINDOW_SIZE_PS} ps ({WINDOW_SIZE_PS / 1000.0} ns)")

    gro_path = os.path.join(folder_path, TOPOLOGY_FILE)
    
    if not os.path.exists(gro_path):
        print(f"Error: Topology file '{TOPOLOGY_FILE}' not found. Exiting.")
        return

    # Find all trajectory files matching the pattern
    xtc_files = sorted([f for f in os.listdir(folder_path) 
                        if f.startswith(TRAJECTORY_PREFIX) and f.endswith(TRAJECTORY_EXT)])
    
    if not xtc_files:
        print(f"No '{TRAJECTORY_PREFIX}*{TRAJECTORY_EXT}' files found in {folder_path}")
        return

    for xtc_file in xtc_files:
        xtc_path = os.path.join(folder_path, xtc_file)
        print(f"\nProcessing trajectory: {xtc_file}")
        
        try:
            u = mda.Universe(gro_path, xtc_path)
            ref = mda.Universe(gro_path)
            n_frames = len(u.trajectory)
            
            if n_frames == 0:
                print("  Error: Trajectory has 0 frames.")
                continue
            print(f"  Trajectory contains {n_frames} frames.")
            
        except Exception as e:
            print(f"  Error loading trajectory {xtc_file}: {e}")
            continue

        # --- STEP 1: Prepare Atom Selections (Sliding Window) ---
        print("  Building atom selections...")
        
        try:
            protein_residues = list(u.select_atoms("protein").residues)
            ref_protein_residues = list(ref.select_atoms("protein").residues)
            
            if not protein_residues:
                print("  Error: No 'protein' selection found.")
                continue
        except Exception as e:
            print(f"  Error selecting protein atoms: {e}")
            continue
            
        traj_selections = []
        ref_selections = []   # Stores reference positions (numpy arrays)
        middle_resid_ids = [] # Maps index to the central Residue ID

        # Create a sliding window of 3 residues
        for i in range(1, len(protein_residues) - 1):
            traj_triple = protein_residues[i-1:i+2]
            ref_triple = ref_protein_residues[i-1:i+2]

            try:
                # Create selection string dynamically
                traj_sel_str = " or ".join([f"resid {res.resid}" for res in traj_triple])
                traj_ag = u.select_atoms(f"backbone and ({traj_sel_str})")
                
                ref_sel_str = " or ".join([f"resid {res.resid}" for res in ref_triple])
                ref_ag = ref.select_atoms(f"backbone and ({ref_sel_str})")

                # Ensure selections match and are valid
                if len(traj_ag) > 0 and len(traj_ag) == len(ref_ag):
                    traj_selections.append(traj_ag)
                    ref_selections.append(ref_ag.positions) 
                    middle_resid_ids.append(protein_residues[i].resid)
                
            except Exception:
                continue
        
        if not traj_selections:
            print("  Error: Could not create any valid triple selections.")
            continue
            
        n_triples = len(traj_selections)
        print(f"  Created {n_triples} valid 3-residue backbone selections.")

        # --- STEP 2: Initialize Result Matrices ---
        all_rmsd_data = np.zeros((n_frames, n_triples), dtype=np.float32)
        frame_times = np.zeros(n_frames, dtype=np.float32)
        
        # --- STEP 3: Trajectory Iteration ---
        print(f"  Starting single pass over {n_frames} frames... (This may take some time)")
        
        for frame_index, ts in enumerate(u.trajectory):
            frame_times[frame_index] = ts.time
            
            if frame_index > 0 and frame_index % 500 == 0:
                print(f"    ...processed frame {frame_index}/{n_frames}", end='\r')

            for triple_index in range(n_triples):
                traj_pos = traj_selections[triple_index].positions
                ref_pos = ref_selections[triple_index]
                
                # Calculate RMSD with superposition
                current_rmsd = rms.rmsd(traj_pos, ref_pos, superposition=True)
                all_rmsd_data[frame_index, triple_index] = current_rmsd
        
        print("\n  ...Trajectory pass completed.")

        # --- STEP 4: Time Binning and Averaging ---
        print("  Calculating binned averages and writing output files...")
        
        total_time = frame_times[-1]
        time_bins = np.arange(0, total_time + WINDOW_SIZE_PS, WINDOW_SIZE_PS)
        
        # Determine filename suffix
        if xtc_file == "prod.xtc":
            traj_suffix_base = "prod"
        else:
            traj_suffix_base = xtc_file.replace(".xtc", "").split("_")[-1]
            
        print(f"  Total time: {total_time:.2f} ps. Creating {len(time_bins) - 1} CSV files.")

        # --- STEP 5: Write CSV Files ---
        for bin_idx in range(len(time_bins) - 1):
            bin_start = time_bins[bin_idx]
            bin_end = time_bins[bin_idx+1]
            
            start_ns = bin_start / 1000.0
            end_ns = bin_end / 1000.0

            # Select frames within the current time window
            mask = (frame_times >= bin_start) & (frame_times < bin_end)
            
            if not np.any(mask):
                print(f"  Warning: No frames found for window {start_ns:.1f}-{end_ns:.1f}ns. Skipping.")
                continue
            
            # Calculate mean RMSD over the time axis
            rmsd_data_in_bin = all_rmsd_data[mask, :]
            avg_rmsd_for_triples = np.mean(rmsd_data_in_bin, axis=0)

            # Define output filename
            csv_filename = f"rmsd_data_{traj_suffix_base}_{start_ns:.1f}-{end_ns:.1f}ns.csv"
            csv_file_path = os.path.join(folder_path, csv_filename)

            try:
                with open(csv_file_path, "w", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["Residue_ID", "Mean_RMSD_nm"])
                    
                    for triple_index, resid in enumerate(middle_resid_ids):
                        rmsd_value = avg_rmsd_for_triples[triple_index]
                        writer.writerow([resid, f"{rmsd_value:.6f}"])
                
                print(f"    Saved: {csv_filename}")
            
            except Exception as e:
                print(f"  Error writing CSV {csv_filename}: {e}")

    print("\nAll processing finished successfully.")

if __name__ == "__main__":
    main()