# python rmsf_graphs.py <pdb_name> <output_path>
# Creates five subplots showing RMSF values for each predictor

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D

def make_rmsf_scatter(pdb_rmsf_list, output_path):
    print(f"[rmsf_graphs.py] Making subplot for {len(pdb_rmsf_list)} PDBs...")
    
    if not pdb_rmsf_list:
        print("[rmsf_graphs.py] No data to plot")
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
    
    for i, pdb_data in enumerate(pdb_rmsf_list):
        pdb_id = pdb_data["pdb_id"]
        rmsf_df = pdb_data["rmsf"]
        bl_rmsf_df = pdb_data["rmsf_bl"]
        
        
        rmsf_df = rmsf_df[rmsf_df['residue'] != 0]
        
        ax.plot(
                    rmsf_df['residue'].drop_duplicates(),
                    bl_rmsf_df['RMSF'],
                    linestyle='-',
                    color="#8f8f8f",
                    alpha=1,
                    lw=2.5,
                    label=f'{pdb_id.upper()} BL' if i == 0 else None,
                    zorder=3
                )
    
        for predictor in predictors.keys():
            predictor_data = rmsf_df[rmsf_df['predictor'] == predictor]
            
            if not predictor_data.empty:
                ax.plot(
                    predictor_data['residue'],
                    predictor_data['rmsf'],
                    linestyle='--',
                    color=predictorcolors.get(predictor),
                    alpha=1,
                    lw=2.5,
                    label=pdb_id.upper() if i == 0 else None
                )

    
    #ax.set_title(f'RMSF by Residue of Ensemble Predictions - {pdb_id.upper()}', fontsize=36)
    ax.set_xlabel('Residue Number', fontsize=42)
    ax.set_ylabel('RMSF (Å)', fontsize=42)
    
    ax.tick_params(axis='both', which='major', labelsize=24)
    
    for _, predictor in enumerate(predictors.keys()):
        pdb_legend_elements.append(
            Line2D([0], [0], color=predictorcolors[predictor], lw=2, label=predictor.capitalize(), linestyle="--")
        )
    pdb_legend_elements.append(
        Line2D([0], [0], color="#8f8f8f", lw=2, label="Deposited", linestyle='-')
    )
    ax.legend(handles=pdb_legend_elements, loc='upper right', bbox_to_anchor=(1.005, 1.02), ncol=len(pdb_rmsf_list),framealpha=0.7, fontsize=28)
    
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"[rmsf_graphs.py] Saved subplot visualization to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsf_graphs.py <pdb_name> <output_path>")
        sys.exit(1)

    pdb_name = sys.argv[1]
    output_path = sys.argv[2]
    pdb_rmsf_set = []

    # modified from dataset to just one pdb
    pdb_list = [ pdb_name ]
    
    for pdb in pdb_list:
        try:
            rmsf_path = f"./PDBs/{pdb}/analysis/rmsf.csv"
            rmsf_bl_path = f"./PDBs/{pdb}/analysis/{pdb}_qfit_RMSF.csv"
            
            if os.path.exists(rmsf_path):
                rmsf_df = pd.read_csv(rmsf_path)
                bl_rmsf_df = pd.read_csv(rmsf_bl_path)
                
                pdb_rmsf_set.append({
                    "pdb_id": pdb,
                    "rmsf": rmsf_df,
                    'rmsf_bl': bl_rmsf_df
                })
            else:
                print(f"[rmsf_graphs.py] RMSF file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsf_graphs.py] Error processing {pdb}: {e}")

    make_rmsf_scatter(pdb_rmsf_set, output_path)