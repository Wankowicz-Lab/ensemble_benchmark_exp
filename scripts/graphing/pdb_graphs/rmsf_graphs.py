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
        
    fig, axes = plt.subplots(5, 1, figsize=(16, 30), sharex=True)

    predictorcolors = {
        'bioemu': '#3498db',    # Blue
        'alphaflow': '#e74c3c', # Red
        'sam2': '#2ecc71',       # Green
        'openfold': '#f1c40f',  # Yellow
        'boltz2': '#9b59b6',  # Purple
    }

    # note: scatter markers are removed for now for testing
    predictors = {
        'sam2': 's',     # Square
        'alphaflow': 'o', # Circle
        'bioemu': '^',   # Triangle
        'openfold': 'D', # Diamond
        'boltz2': 'X'    # X-mark
    }
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(pdb_rmsf_list)))
    pdb_legend_elements = []
    
    predictor_to_ax = {
        'sam2': axes[0],
        'alphaflow': axes[1],
        'bioemu': axes[2],
        'openfold': axes[3],
        'boltz2': axes[4]
    }
    
    # Set titles for each subplot
    for predictor, ax in predictor_to_ax.items():
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_ylabel('RMSF (Å)', fontsize=24)
        ax.set_ylim(bottom=0)
        
        if ax != axes[-1]:
            ax.tick_params(axis='x', labelbottom=False)
    
    # Only set x-label on the bottom subplot
    axes[-1].set_xlabel('Residue', fontsize=24)
    
    for i, pdb_data in enumerate(pdb_rmsf_list):
        pdb_id = pdb_data["pdb_id"]
        rmsf_df = pdb_data["rmsf"]
        color = colors[i]
        
        for predictor in predictors.keys():
            predictor_data = rmsf_df[rmsf_df['predictor'] == predictor]
            
            if not predictor_data.empty:
                ax = predictor_to_ax[predictor]
                
                ax.plot(
                    predictor_data['residue'],
                    predictor_data['rmsf'],
                    linestyle='-',
                    color=predictorcolors.get(predictor),
                    alpha=1,
                    lw=2.5,
                    label=pdb_id.upper() if i == 0 else None
                )
    
    fig.suptitle(f'RMSF by Residue of Ensemble Predictions - {pdb_id.upper()}', fontsize=30)
    
    for _, predictor in enumerate(predictors.keys()):
        pdb_legend_elements.append(
            Line2D([0], [0], color=predictorcolors[predictor], lw=2, label=predictor.capitalize())
        )

    fig.legend(handles=pdb_legend_elements, loc='upper center', 
               bbox_to_anchor=(0.915, 0.085), ncol=len(pdb_rmsf_list),
               title='Predictors', framealpha=0.7, fontsize=20)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.1, hspace=0.05) 
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
            
            if os.path.exists(rmsf_path):
                rmsf_df = pd.read_csv(rmsf_path)
                
                pdb_rmsf_set.append({
                    "pdb_id": pdb,
                    "rmsf": rmsf_df
                })
            else:
                print(f"[rmsf_graphs.py] RMSF file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsf_graphs.py] Error processing {pdb}: {e}")

    make_rmsf_scatter(pdb_rmsf_set, output_path)