# python ediam_graphs.py <pdb_name> <output_path>
# Creates five subplots showing EDIAm values for each predictor

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D
import json

def make_ediam_scatter(pdb_ediam_list, output_path):
    print(f"[ediam_graphs.py] Making subplot for {len(pdb_ediam_list)} PDBs...")
    
    if not pdb_ediam_list:
        print("[ediam_graphs.py] No data to plot")
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
        'openfold': 'D',  # Diamond
        'boltz2': 'X',   # X
    }
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(pdb_ediam_list)))
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
        ax.set_ylabel('EDIAm', fontsize=24)
        ax.set_ylim(bottom=0)
        
        if ax != axes[-1]:
            ax.tick_params(axis='x', labelbottom=False)
    
    # Only set x-label on the bottom subplot
    axes[-1].set_xlabel('Residue', fontsize=24)

    for i, pdb_data in enumerate(pdb_ediam_list):
        pdb_id = pdb_data["pdb_id"]
        ediam_df = pdb_data["ediam"]
        color = colors[i]
        
        for predictor in predictors.keys():
            predictor_data = ediam_df[ediam_df['predictor'] == predictor]

            if not predictor_data.empty:
                ax = predictor_to_ax[predictor]
                
                ax.plot(
                    predictor_data['residue'],
                    predictor_data['EDIAm'], 
                    linestyle='-',
                    color=predictorcolors.get(predictor),
                    alpha=1,
                    lw=2.5,
                    label=pdb_id.upper() if i == 0 else None
                )

    fig.suptitle(f'EDIAm by Residue of Ensemble Predictions - {pdb_id.upper()}', fontsize=30)
    
    for _, predictor in enumerate(predictors.keys()):
        pdb_legend_elements.append(
            Line2D([0], [0], color=predictorcolors[predictor], lw=2, label=predictor.capitalize())
        )

    fig.legend(handles=pdb_legend_elements, loc='upper center', 
               bbox_to_anchor=(0.915, 0.085), ncol=len(pdb_ediam_list),
               title='Predictors', framealpha=0.7, fontsize=20)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.95, bottom=0.1, hspace=0.05) 
    plt.savefig(output_path)
    print(f"[ediam_graphs.py] Saved subplot visualization to {output_path}")

def parse_density_fitness_csv(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.csv"
    
    if not os.path.exists(file_path):
        print(f"[ediam_graphs.py] File not found: {file_path}") 
        return pd.DataFrame()
    
    try:
        lines = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    lines.append({
                        'predictor': line.split(',')[0].strip(),
                        'frame': line.split(',')[1].strip(),
                        'metrics': (','.join(line.split(',')[2:]).strip())[1:-1]  # Remove surrounding quotes
                    })
        
        df = pd.DataFrame(lines)
        if 'metrics' in df.columns:
            all_data = []
            
            for _, row in df.iterrows():
                predictor = row['predictor']
                frame = row['frame']
                
                try:
                    metrics = json.loads(row['metrics'])
                    
                    for residue_data in metrics:
                        if 'EDIAm' in residue_data:
                            all_data.append({
                                'predictor': predictor,
                                'frame': frame,
                                'residue': residue_data.get('seqID', None),
                                'aa': residue_data.get('compID', None),
                                'RSCCS': residue_data['RSCCS'],
                                'RSR': residue_data.get('RSR', None),
                                'EDIAm': residue_data.get('EDIAm', None),
                                'chain': residue_data.get('asymID', None)
                            })
                except json.JSONDecodeError:
                    print(row['metrics'])
                    print(f"[ediam_graphs.py] Error parsing metrics JSON for {predictor} frame {frame}")
                except Exception as e:
                    print(f"[ediam_graphs.py] Error processing metrics: {e}")
            
            return pd.DataFrame(all_data)
            
    except Exception as e:
        print(f"[ediam_graphs.py] Error reading {file_path}: {e}")
        return pd.DataFrame()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python ediam_graphs.py <pdb_name> <output_path>")
        sys.exit(1)

    pdb_name = sys.argv[1]
    output_path = sys.argv[2]
    pdb_ediam_set = []

    # modified from dataset to just one pdb
    pdb_list = [pdb_name]
    
    for pdb in pdb_list:
        try:
            ediam_df = parse_density_fitness_csv(pdb)
            
            if not ediam_df.empty:
                pdb_ediam_set.append({
                    "pdb_id": pdb,
                    "ediam": ediam_df 
                })
            else:
                print(f"[ediam_graphs.py] EDIAm data not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[ediam_graphs.py] Error processing {pdb}: {e}")

    make_ediam_scatter(pdb_ediam_set, output_path) 