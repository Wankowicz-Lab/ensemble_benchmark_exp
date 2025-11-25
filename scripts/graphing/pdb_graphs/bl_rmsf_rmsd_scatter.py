# python bl_rmsf_rmsd_scatter.py <pdb_name> <output_path>
# Creates a graph showing baseline PDB QFIT RMSF as a line plot with predictor RMSF as colored scatter points

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D


def read_rmsr_values(pdb_id):
    pdb_dir = "./PDBs/" + pdb_id
    rmsf_file = pdb_dir + '/analysis/rmsr_galign_each.csv'
    
    if not os.path.exists(rmsf_file):
        print(f"[bl_rmsf_rmsd_scatter.py] Warning: {rmsf_file} not found")
        return pd.DataFrame()
    
    thedf = pd.read_csv(rmsf_file)
    return thedf


def get_deposited_rmsf_values(pdb_id):
    pdb_dir = "./PDBs/" + pdb_id
    rmsf_file = pdb_dir + '/analysis/' + pdb_id + "_qfit_RMSF.csv"
    
    if not os.path.exists(rmsf_file):
        print(f"[bl_rmsf_rmsd_scatter.py] Warning: {rmsf_file} not found")
        return pd.DataFrame()
    
    thedf = pd.read_csv(rmsf_file)
    return thedf


def create_baseline_rmsf_scatter(pdb_id, rmsr_df, deposited_rmsf, output_path):
    
    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c',  # Red
        'bioemu': '#3498db',     # Blue
        'openfold': '#f1c40f',   # Yellow
        'boltz2': '#9b59b6',     # Purple
    }
    
    predictor_markers = {
        'sam2': 'o',      # Circle
        'alphaflow': 's',  # Square
        'bioemu': '^',     # Triangle up
        'openfold': 'D',   # Diamond
        'boltz2': 'v',     # Triangle down
    }
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    rmsr_df = rmsr_df[rmsr_df['residue'] != 0]
    
    if 'RMSF' in deposited_rmsf.columns:
        baseline_rmsf_col = 'RMSF'
    elif 'rmsf' in deposited_rmsf.columns:
        baseline_rmsf_col = 'rmsf'
    else:
        print(f"[bl_rmsf_rmsd_scatter.py] Error: Could not find RMSF column in baseline data")
        return
    
    baseline_residues = sorted(deposited_rmsf['residue'].unique()) if 'residue' in deposited_rmsf.columns else sorted(deposited_rmsf.index)
    baseline_rmsf_values = deposited_rmsf[baseline_rmsf_col].values
    
    ax.plot(baseline_residues, baseline_rmsf_values, 
            linestyle='-', color='#2c3e50', linewidth=3, 
            label='Baseline Multiconformer RMSF', zorder=5, alpha=0.8)
    
    for predictor in predictorcolors.keys():
        predictor_data = rmsr_df[rmsr_df['predictor'] == predictor]
        
        if not predictor_data.empty:
            ax.scatter(predictor_data['residue'], predictor_data['RMSR'],
                      color=predictorcolors[predictor], 
                      marker='x',
                      s=60, alpha=1, edgecolors='black', linewidth=0.5,
                      label=f'{predictor.upper()}', zorder=4)
    
    ax.set_xlabel('Residue Number', fontsize=28)
    ax.set_ylabel('RMSF/RMSR (Å)', fontsize=28)

    
    ax.grid(True, linestyle='--', alpha=0.3, zorder=1)
    
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    ax.set_ylim(bottom=0)
    
    legend_elements = [Line2D([0], [0], color='#2c3e50', lw=3, label='Baseline Multiconformer RMSF')]
    
    for predictor in predictorcolors.keys():
        predictor_data = rmsr_df[rmsr_df['predictor'] == predictor]
        if not predictor_data.empty:
            legend_elements.append(
                Line2D([0], [0], marker='o', 
                      color='w', markerfacecolor=predictorcolors[predictor],
                      markersize=8, markeredgecolor='black', markeredgewidth=0.5,
                      label=f'{predictor.upper()}', linestyle='None')
            )
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=16, 
              frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[bl_rmsf_rmsd_scatter.py] Saved plot to {output_path}")
    plt.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python bl_rmsf_rmsd_scatter.py <pdb_name> <output_path>")
        sys.exit(1)

    pdb_name = sys.argv[1]
    output_path = sys.argv[2]

    rmsf_file = read_rmsr_values(pdb_name)
    deposited_rmsf = get_deposited_rmsf_values(pdb_name)
    
    if rmsf_file.empty:
        print(f"[bl_rmsf_rmsd_scatter.py] Error: No predictor RMSF data found for {pdb_name}")
        sys.exit(1)
    
    if deposited_rmsf.empty:
        print(f"[bl_rmsf_rmsd_scatter.py] Error: No baseline QFIT RMSF data found for {pdb_name}")
        sys.exit(1)

    create_baseline_rmsf_scatter(pdb_name, rmsf_file, deposited_rmsf, output_path)
    
    print(f"[bl_rmsf_rmsd_scatter.py] Successfully created baseline RMSF scatter plot for {pdb_name}")