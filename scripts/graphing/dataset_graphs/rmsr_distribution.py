# python rmsr_distribution.py <dataset_name> <output_path>
# Makes a graph of the rmsr distribution for a given dataset with all PDBs.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D


def make_dataset_rmsr_boxplot(pdb_rmsr_list, output_path, dataset_name):
    print(f"[rmsr_distribution.py] Making rmsr box plot for {len(pdb_rmsr_list)} PDBs...")
    
    if not pdb_rmsr_list:
        print("[rmsr_distribution.py] No data to plot")
        return

    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c', # Red
        'bioemu': '#3498db',    # Blue
        'openfold': '#f1c40f',  # Yellow
        'boltz2': '#9b59b6',  # Purple
        'Deposited': '#95a5a6'    # Gray
    }

    predictor_order = ['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2']
    
    box_data = []
    box_labels = []
    box_colors = []
    positions = []
    group_boundaries = []  # Track where to place vertical lines
    
    current_pos = 1
    group_spacing = 1.2
    box_spacing = 0.8
    
    for i, pdb_data in enumerate(pdb_rmsr_list):
        pdb_id = pdb_data["pdb_id"]
        rmsr_df = pdb_data["rmsr"]
        deposited_rmsr = pdb_data.get("rmsr_bl")
        

        predictor_rmsrs = {
            'sam2': rmsr_df[rmsr_df['predictor'] == 'sam2'],
            'alphaflow': rmsr_df[rmsr_df['predictor'] == 'alphaflow'],
            'bioemu': rmsr_df[rmsr_df['predictor'] == 'bioemu'],
            'openfold': rmsr_df[rmsr_df['predictor'] == 'openfold'],
            'boltz2': rmsr_df[rmsr_df['predictor'] == 'boltz2'],
        }

        if deposited_rmsr is not None:
            deposited_rmsr['rmsr'] = deposited_rmsr['RMSR']
            predictor_rmsrs['Deposited'] = deposited_rmsr[deposited_rmsr['rmsr'] > 0]


        for j, predictor in enumerate(predictor_order):
            data = predictor_rmsrs.get(predictor)
            if data is not None and not data.empty:
                box_data.append(data['RMSR'].values)
                box_colors.append(predictorcolors[predictor])
                positions.append(current_pos + j * box_spacing)
                
                if i == 0:
                    box_labels.append(predictor.upper())
                else:
                    box_labels.append("")



        if i < len(pdb_rmsr_list) - 1:
            group_end = current_pos + (len(predictor_order) - 1) * box_spacing + group_spacing / 2
            group_boundaries.append(group_end + box_spacing / 2)
        
        current_pos += len(predictor_order) * box_spacing + group_spacing

    if not box_data:
        print("[rmsr_distribution.py] No valid data found")
        return

    fig, ax = plt.subplots(figsize=(24, 10))
    
    bp = ax.boxplot(box_data, positions=positions, patch_artist=True, 
                    showmeans=True, meanline=True, showfliers=False,
                    widths=0.6)
     
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
    
    plt.setp(bp['means'], color='black', linewidth=1)
    
    for boundary in group_boundaries:
        ax.axvline(x=boundary, color='gray', linestyle='-', alpha=0.5, linewidth=1)
    
    pdb_centers = []
    current_center = 1 + (len(predictor_order) - 1) * box_spacing / 2
    
    for i, pdb_data in enumerate(pdb_rmsr_list):
        pdb_centers.append(current_center)
        current_center += len(predictor_order) * box_spacing + group_spacing
    
    ax.set_xticks(pdb_centers)
    ax.set_xticklabels([pdb_data["pdb_id"].upper() for pdb_data in pdb_rmsr_list], 
                        fontsize=24)
    
    ax.set_yticklabels(ax.get_yticks(), fontsize=24)
    
    legend_elements = []
    for predictor in predictor_order:
        legend_elements.append(
            plt.Rectangle((0,0),1,1, facecolor=predictorcolors[predictor], 
                         alpha=0.7, label=predictor.upper())
        )
    

    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.1), 
              fontsize=28, framealpha=0.9)
    
    ax.set_ylabel('RMSR (Å)', fontsize=36)
    ax.set_xlabel('PDB Structure', fontsize=36)
    # ax.set_title(f'rmsr Distribution by Predictor - {dataset_name.upper()} Dataset', 
    #              fontsize=18)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"[rmsr_distribution.py] rmsr distribution box plot saved to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsr_distribution.py <dataset_name> <output_path>")
        sys.exit(1)

    dataset_name = sys.argv[1]
    output_path = sys.argv[2]
    
    pdb_rmsr_set = []

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
            rmsr_path = f"./PDBs/{pdb}/analysis/rmsr_galign.csv"
            
            if os.path.exists(rmsr_path):
                rmsr_df = pd.read_csv(rmsr_path)
                
                pdb_data = {
                    "pdb_id": pdb,
                    "rmsr": rmsr_df
                }
                
                pdb_rmsr_set.append(pdb_data)
            else:
                print(f"[rmsr_distribution.py] rmsr file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsr_distribution.py] Error processing {pdb}: {e}")

    make_dataset_rmsr_boxplot(pdb_rmsr_set, output_path, dataset_name)