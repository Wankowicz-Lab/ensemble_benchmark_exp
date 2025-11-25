# python rmsr_condensed.py <pdb_name> <output_path>
# uses global aligned data analysis/rmsr_galign.csv from get_rmsr_galign.py, which aligned with pymol.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D

def make_rmsr_scatter(pdb_rmsr_list, output_path):
    print(f"[rmsr_condensed.py] Making subplot for {len(pdb_rmsr_list)} PDBs...")
    
    if not pdb_rmsr_list:
        print("[rmsr_condensed.py] No data to plot")
        return
        
    predictorcolors = {
        'bioemu': '#3498db',    # Blue
        'alphaflow': '#e74c3c', # Red
        'sam2': '#2ecc71',       # Green
        'openfold': '#f1c40f',  # Yellow
        'boltz2': '#9b59b6',  # Purple
    }

    predictors = {
        'sam2': 's',     # Square
        'alphaflow': 'o', # Circle
        'bioemu': '^',   # Triangle
        'openfold': 'D', # Diamond
        'boltz2': 'X'    # X-mark
    }
    
    pdb_legend_elements = []
    

    fig, ax = plt.subplots(figsize=(24, 8))
    
    for i, pdb_data in enumerate(pdb_rmsr_list):
        pdb_id = pdb_data["pdb_id"]
        rmsr_df = pdb_data["rmsr"]
        
        for predictor in predictors.keys():
            predictor_data = rmsr_df[rmsr_df['predictor'] == predictor]
            
            if not predictor_data.empty:
                ax.plot(
                    predictor_data['residue'],
                    predictor_data['RMSR'],
                    linestyle='-',
                    color=predictorcolors.get(predictor),
                    alpha=1,
                    lw=2.5,
                    label=pdb_id.upper() if i == 0 else None
                )

    
    ax.set_title(f'RMSR by Residue of Ensemble Predictions from Deposited PDB - {pdb_id.upper()}', fontsize=36)
    ax.set_xlabel('Residue Number', fontsize=42)
    ax.set_ylabel('RMSR (Å)', fontsize=42)
    
    ax.tick_params(axis='both', which='major', labelsize=24)
    
    for _, predictor in enumerate(predictors.keys()):
        pdb_legend_elements.append(
            Line2D([0], [0], color=predictorcolors[predictor], lw=2, label=predictor.capitalize(), linestyle="-")
        )
    ax.legend(handles=pdb_legend_elements, loc='upper right', bbox_to_anchor=(1.005, 1.02), ncol=len(pdb_rmsr_list),framealpha=0.7, fontsize=28)
    
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"[rmsr_condensed.py] Saved subplot visualization to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsr_condensed.py <pdb_name> <output_path>")
        sys.exit(1)

    pdb_name = sys.argv[1]
    output_path = sys.argv[2]
    pdb_rmsr_set = []

    # modified from dataset to just one pdb
    pdb_list = [ pdb_name ]
    
    for pdb in pdb_list:
        try:
            rmsr_path = f"./PDBs/{pdb}/analysis/rmsr_galign.csv"
            
            if os.path.exists(rmsr_path):
                rmsr_df = pd.read_csv(rmsr_path)
                
                pdb_rmsr_set.append({
                    "pdb_id": pdb,
                    "rmsr": rmsr_df
                })
            else:
                print(f"[rmsr_condensed.py] rmsr file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsr_condensed.py] Error processing {pdb}: {e}")

    make_rmsr_scatter(pdb_rmsr_set, output_path)