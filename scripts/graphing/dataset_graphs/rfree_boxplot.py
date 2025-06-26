# python dataset_rfree_boxplot.py <dataset_name> <output_path>
# Makes a box plot of R-Free values from the predictors

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


def make_boxplot(pdb_rfrees_list, output_path):
    print(f"[dataset_rfree_boxplot.py] Making box plot for {len(pdb_rfrees_list)} PDBs...")
    
    if not pdb_rfrees_list:
        print("[dataset_rfree_boxplot.py] No data to plot")
        return
    
    predictor_data = {
        'sam2': [],
        'alphaflow': [],
        'bioemu': [],
        'openfold': [],
        'boltz2': [],
    }
    
    for pdb_data in pdb_rfrees_list:
        rfrees_df = pdb_data["rfrees"]
        
        for predictor in predictor_data.keys():
            predictor_rows = rfrees_df[rfrees_df['predictor'] == predictor]
            if not predictor_rows.empty:
                predictor_data[predictor].append(predictor_rows['rfree'].iloc[0])
    
    data_to_plot = []
    labels = []
    colors = ['#2ecc71', '#e74c3c', '#3498db', '#f1c40f', '#9b59b6' ]  # sam2, alphaflow, bioemu, openfold. boltz2

    for predictor in ['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2']:
        if predictor_data[predictor]:
            data_to_plot.append(predictor_data[predictor])
            labels.append(predictor.upper())
    
    if not data_to_plot:
        print("[dataset_rfree_boxplot.py] No valid data found for any predictor")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    boxprops = {'alpha': 0.7}
    medianprops = {'color': 'black', 'linestyle': '--', 'linewidth': 1}
    meanprops = {'color': (32/255,32/255,32/255), 'linestyle': '--', 'linewidth': 1}

    box_plot = ax.boxplot(data_to_plot, 
                         labels=labels,
                         patch_artist=True,
                         showmeans=True,
                         meanline=True,
                         meanprops=meanprops,
                         boxprops=boxprops,
                         medianprops=medianprops)
    
    for patch, color in zip(box_plot['boxes'], colors[:len(data_to_plot)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_ylabel('R-free', fontsize=14)
    ax.set_xlabel('Predictor', fontsize=14)
    ax.set_title('R-Free Distribution Comparison of Ensemble Predictions', fontsize=16)
    
    # Set y-axis limits (R-free is between 0 and 1)
    ylim_bottom = max(0, ax.get_ylim()[0] - 0.05)
    ylim_top = min(1.0, ax.get_ylim()[1] + 0.05)
    ax.set_ylim(ylim_bottom, ylim_top)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"[dataset_rfree_boxplot.py] Saved box plot to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dataset_rfree_boxplot.py <dataset_name> <output_path>")
        sys.exit(1)

    dataset_name = sys.argv[1]
    output_path = sys.argv[2]
    pdb_rfrees_set = []

    pdb_list = []
    with open(f"./splits/{dataset_name}.txt") as f:
        for line in f:
            pdb_list.append(line.strip())
    
    for pdb in pdb_list:
        try:
            rfree_path = f"./PDBs/{pdb}/analysis/rfrees.csv"
            rfree_df = pd.read_csv(rfree_path)
            
            pdb_rfrees_set.append({
                "pdb_id": pdb,
                "rfrees": rfree_df
            })
            
        except FileNotFoundError:
            print(f"[dataset_rfree_boxplot.py] File not found for {pdb}, skipping...")

    make_boxplot(pdb_rfrees_set, output_path)