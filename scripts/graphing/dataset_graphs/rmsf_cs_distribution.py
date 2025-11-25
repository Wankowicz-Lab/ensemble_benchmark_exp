# python rmsf_cs_distribution.py <dataset_name> <output_path>
# Makes a graph of the rmsf cosine similarity distribution for a given dataset with all PDBs.
# the "all" column is for deposited only, not predictor to predictor comparison

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D


def make_cosine_similarity_boxplot(pdb_cs_list, output_path, dataset_name):
    print(f"[rmsf_cs_distribution.py] Making cosine similarity box plot for {len(pdb_cs_list)} PDBs...")
    
    if not pdb_cs_list:
        print("[rmsf_cs_distribution.py] No data to plot")
        return

    # Define colors for the predictors
    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c',   # Red
        'bioemu': '#3498db',      # Blue
        'openfold': '#f1c40f',    # Yellow
        'boltz2': '#9b59b6',      # Purple
        'All': '#95a5a6'          # Gray
    }

    predictor_order = ['All', 'sam2', 'alphaflow', 'bioemu', 'boltz2', 'openfold']
    
    all_deposited_cs_values = []
    predictor_cs_values = {
        'sam2': [],
        'alphaflow': [],
        'bioemu': [],
        'boltz2': [],
        'openfold': []
    }
    
    for pdb_data in pdb_cs_list:
        cs_df = pdb_data["cosine_similarity"]
        
        deposited_comparisons = cs_df[(cs_df['element1'] == 'deposited') | (cs_df['element2'] == 'deposited')]
        
        all_deposited_cs_values.extend(deposited_comparisons['cosine_similarity'].tolist())
        
        for predictor in predictor_cs_values.keys():
            predictor_rows = deposited_comparisons[
                ((deposited_comparisons['element1'] == predictor) & (deposited_comparisons['element2'] == 'deposited')) |
                ((deposited_comparisons['element1'] == 'deposited') & (deposited_comparisons['element2'] == predictor))
            ]
            if not predictor_rows.empty:
                predictor_cs_values[predictor].extend(predictor_rows['cosine_similarity'].tolist())
    
    box_data = [all_deposited_cs_values] + [predictor_cs_values[pred] for pred in predictor_order if pred != 'All']
    box_colors = [predictorcolors[pred] for pred in predictor_order]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    bp = ax.boxplot(box_data, patch_artist=True, 
                    showmeans=True, meanline=True, showfliers=False,
                    widths=0.6)
    
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
    
    plt.setp(bp['means'], color='black', linewidth=1)
    
    ax.set_xticklabels([pred.upper() for pred in predictor_order], fontsize=16)
    ax.set_yticklabels([f"{tick:.2f}" for tick in ax.get_yticks()], fontsize=16)
    
    ax.set_ylabel('Cosine Similarity', fontsize=24)
    #ax.set_title(f'Cosine Similarity Distribution - {dataset_name}', fontsize=24)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    print(f"[rmsf_cs_distribution.py] Cosine similarity distribution plot saved to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsf_distribution.py <dataset_name> <output_path>")
        sys.exit(1)

    dataset_name = sys.argv[1]
    output_path = sys.argv[2]

    pdb_rmsf_set = []

    pdb_list = []
    try:
        with open(f"./splits/{dataset_name}.txt") as f:
            for line in f:
                pdb_list.append(line.strip())
    except FileNotFoundError:
        print(f"Dataset file ./splits/{dataset_name}.txt not found")
        sys.exit(1)
    
    for pdb in pdb_list:
        try:
            cosine_similarity_path = f"./PDBs/{pdb}/analysis/cosine_similarity.csv"
            if os.path.exists(cosine_similarity_path):
                cosine_similarity = pd.read_csv(cosine_similarity_path)
                
                pdb_data = {
                    "pdb_id": pdb,
                    "cosine_similarity": cosine_similarity
                }
                
                pdb_rmsf_set.append(pdb_data)
            else:
                print(f"[rmsf_cs_distribution.py] RMSF file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsf_cs_distribution.py] Error processing {pdb}: {e}")

    make_cosine_similarity_boxplot(pdb_rmsf_set, output_path, dataset_name)