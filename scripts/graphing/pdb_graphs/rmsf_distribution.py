# python rmsf_distribution.py <pdb_name> <output_path>
# Makes a graph of the rmsf distribution for a given PDB file.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D


def read_rmsf_values(pdb_id):
    pdb_dir = "./PDBs/" + pdb_id
    rmsf_file = pdb_dir + '/analysis/rmsf.csv'
    thedf = pd.read_csv(rmsf_file)
    return thedf

def get_deposited_rmsf_values(pdb_id):
    pdb_dir = "./PDBs/" + pdb_id
    rmsf_file = pdb_dir + '/analysis/' + pdb_id + "_qfit_RMSF.csv"
    thedf = pd.read_csv(rmsf_file)
    return thedf

if (__name__ == "__main__"):
    if len(sys.argv) != 3:
        print("Usage: python rmsf_distribution.py <pdb_name> <output_path>")
        sys.exit(1)

    pdb_name = sys.argv[1]
    output_path = sys.argv[2]

    rmsf_file = read_rmsf_values(pdb_name)
    deposited_rmsf = get_deposited_rmsf_values(pdb_name)
    deposited_rmsf['rmsf'] = deposited_rmsf['RMSF']


    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c', # Red
        'bioemu': '#3498db',    # Blue
        'openfold': '#f1c40f',  # Yellow
        'boltz2': '#9b59b6',  # Purple
        'Deposited': '#95a5a6'    # Gray
    }


    predictor_rmsfs = {
        'sam2': rmsf_file[rmsf_file['predictor'] == 'sam2'],
        'alphaflow': rmsf_file[rmsf_file['predictor'] == 'alphaflow'],
        'bioemu': rmsf_file[rmsf_file['predictor'] == 'bioemu'],
        'openfold': rmsf_file[rmsf_file['predictor'] == 'openfold'],
        'boltz2': rmsf_file[rmsf_file['predictor'] == 'boltz2'],
        'Deposited': deposited_rmsf[ deposited_rmsf['rmsf'] > 0 ]
    }

    fig, ax = plt.subplots(figsize=(12, 8))

    box_data = []
    box_labels = []
    box_colors = []
    
    for predictor, data in predictor_rmsfs.items():
        if not data.empty:
            box_data.append(data['rmsf'].values)
            box_labels.append(predictor.upper())
            box_colors.append(predictorcolors[predictor])
    
    bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, 
                    showmeans=True, meanline=True, showfliers=False)
    
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
    
    plt.setp(bp['means'], color='black', linewidth=1)
    
    ax.set_xlabel('Predictor', fontsize=14)
    ax.set_ylabel('RMSF (Å)', fontsize=14)
    ax.set_title(f'RMSF Distribution by Predictor - {pdb_name.upper()}', fontsize=16)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"RMSF distribution box plot saved to: {output_path}")




