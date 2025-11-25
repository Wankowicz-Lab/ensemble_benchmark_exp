# python rmsf_distribution.py <dataset_name> <output_path>
# Makes a graph of the rmsf distribution for a given dataset with all PDBs.

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

def read_protein_class_mapping(mapping_text=None, mapping_file="bin/protein_classes/include.txt"):
    mapping = {}
    
    if mapping_text:
        lines = mapping_text.strip().split('\n')
        for line in lines:
            if line.strip():
                protein_class, pdb = line.strip().split(',')
                if protein_class not in mapping:
                    mapping[protein_class] = []
                mapping[protein_class].append(pdb)
    elif os.path.exists(mapping_file):
        try:
            with open(mapping_file, 'r') as f:
                for line in f:
                    if line.strip():
                        protein_class, pdb = line.strip().split(',')
                        if protein_class not in mapping:
                            mapping[protein_class] = []
                        mapping[protein_class].append(pdb)
        except Exception as e:
            print(f"Error reading protein class mapping: {e}")
    
    return mapping


def make_dataset_rmsf_boxplot(pdb_rmsf_list, output_path, dataset_name, protein_class_mapping=None):
    print(f"[rmsf_distribution.py] Making RMSF box plot for {len(pdb_rmsf_list)} PDBs...")
    
    if not pdb_rmsf_list:
        print("[rmsf_distribution.py] No data to plot")
        return

    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c',   # Red
        'bioemu': '#3498db',      # Blue
        'openfold': '#f1c40f',    # Yellow
        'boltz2': '#9b59b6',      # Purple
        'Deposited': '#95a5a6'    # Gray
    }

    predictor_order = ['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2', 'Deposited']
    
    pdb_ignore_list = ['9ewk']

    class_grouped_pdbs = {}
    
    if protein_class_mapping:
        pdb_to_class = {}
        for cls, pdbs in protein_class_mapping.items():
            for pdb in pdbs:
                pdb_to_class[pdb.lower()] = cls
        
        for pdb_data in pdb_rmsf_list:
            pdb_id = pdb_data["pdb_id"].lower()

            if pdb_id in pdb_ignore_list:
                continue

            protein_class = pdb_to_class.get(pdb_id, "Other")
            
            if protein_class not in class_grouped_pdbs:
                class_grouped_pdbs[protein_class] = []
            class_grouped_pdbs[protein_class].append(pdb_data)
    else:
        class_grouped_pdbs["Other"] = pdb_rmsf_list
    
    box_data = []
    box_labels = []
    box_colors = []
    positions = []
    group_boundaries = []
    class_boundaries = []
    class_centers = []
    class_names = []
    
    current_pos = 1
    group_spacing = 1.2
    box_spacing = 0.8
    class_spacing = 3
    
    print (class_grouped_pdbs.keys())
    sorted_classes = ["Lysozyme", "CypA", "UBI", "STN", "Other"]#(class_grouped_pdbs.keys())
    
    for class_idx, protein_class in enumerate(sorted_classes):
        pdbs_in_class = class_grouped_pdbs[protein_class]
        
        class_start_pos = current_pos
        
        for i, pdb_data in enumerate(pdbs_in_class):
            pdb_id = pdb_data["pdb_id"]
            rmsf_df = pdb_data["rmsf"]
            deposited_rmsf = pdb_data.get("rmsf_bl")
            
            predictor_rmsfs = {
                'sam2': rmsf_df[rmsf_df['predictor'] == 'sam2'],
                'alphaflow': rmsf_df[rmsf_df['predictor'] == 'alphaflow'],
                'bioemu': rmsf_df[rmsf_df['predictor'] == 'bioemu'],
                'openfold': rmsf_df[rmsf_df['predictor'] == 'openfold'],
                'boltz2': rmsf_df[rmsf_df['predictor'] == 'boltz2'],
            }
            
            if deposited_rmsf is not None:
                deposited_rmsf['rmsf'] = deposited_rmsf['RMSF']
                predictor_rmsfs['Deposited'] = deposited_rmsf[deposited_rmsf['rmsf'] > 0]
            
            for j, predictor in enumerate(predictor_order):
                data = predictor_rmsfs.get(predictor)
                if data is not None and not data.empty:
                    box_data.append(data['rmsf'].values)
                    box_colors.append(predictorcolors[predictor])
                    positions.append(current_pos + j * box_spacing)
                    
                    if class_idx == 0 and i == 0:
                        box_labels.append(predictor.upper())
                    else:
                        box_labels.append("")
            
            if i < len(pdbs_in_class) - 1:
                group_end = current_pos + (len(predictor_order) - 1) * box_spacing + group_spacing / 2
                group_boundaries.append(group_end + box_spacing / 2)
            
            current_pos += len(predictor_order) * box_spacing + group_spacing
        
        class_centers.append((class_start_pos + current_pos - group_spacing) / 2)
        class_names.append(protein_class)
        
        if class_idx < len(sorted_classes) - 1:
            class_boundaries.append(current_pos )
            current_pos += class_spacing

    if not box_data:
        print("[rmsf_distribution.py] No valid data found")
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
        ax.axvline(x=boundary, color='gray', linestyle='-', alpha=0.3, linewidth=1)
    
    for boundary in class_boundaries:
        ax.axvline(x=boundary, color='black', linestyle='--', alpha=0.5, linewidth=2)
    
    pdb_centers = []
    current_center = 1 + (len(predictor_order) - 1) * box_spacing / 2
    
    for class_idx, protein_class in enumerate(sorted_classes):
        pdbs_in_class = class_grouped_pdbs[protein_class]
        
        for i, pdb_data in enumerate(pdbs_in_class):
            pdb_centers.append(current_center)
            current_center += len(predictor_order) * box_spacing + group_spacing
        
        if class_idx < len(sorted_classes) - 1:
            current_center += class_spacing
    
    ax.set_xticks(pdb_centers)
    pdb_labels = []
    for class_idx, protein_class in enumerate(sorted_classes):
        for pdb_data in class_grouped_pdbs[protein_class]:
            pdb_labels.append(pdb_data["pdb_id"].upper())
    ax.set_xticklabels(pdb_labels, fontsize=16)
    
    for center, class_name in zip(class_centers, class_names):
        ax.annotate(class_name, xy=(center, -0.08), xycoords=('data', 'axes fraction'),
                   ha='center', va='center', fontsize=24, fontweight='bold')
    
    ax.set_yticklabels(ax.get_yticks(), fontsize=24)
    
    legend_elements = []
    for predictor in predictor_order:
        legend_elements.append(
            plt.Rectangle((0,0),1,1, facecolor=predictorcolors[predictor], 
                         alpha=0.7, label=predictor.upper())
        )
    
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.1), 
              fontsize=20, framealpha=0.9)
    
    ax.set_ylabel('RMSF (Å)', fontsize=36)
    ax.set_xlabel('', fontsize=36)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    print(f"[rmsf_distribution.py] RMSF distribution box plot saved to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsf_distribution.py <dataset_name> <output_path>")
        sys.exit(1)

    dataset_name = sys.argv[1]
    output_path = sys.argv[2]
    
    mapping_text = """CypA,4yuo
CypA,3k0m
CypA,5f66
CypA,3k0n
CypA,3k0o
CypA,4yuh
CypA,4yuj
Lysozyme,2lzt
Lysozyme,6o2h
Lysozyme,4lzt
Lysozyme,2vb1
Lysozyme,1aki
Lysozyme,8dyz
Lysozyme,8dz7
UBI,5tof
UBI,5tog
STN,2pyk
STN,2pzw"""
    
    protein_class_mapping = read_protein_class_mapping(mapping_text)
    
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
            rmsf_path = f"./PDBs/{pdb}/analysis/rmsf.csv"
            rmsf_bl_path = f"./PDBs/{pdb}/analysis/{pdb}_qfit_RMSF.csv"
            
            if os.path.exists(rmsf_path):
                rmsf_df = pd.read_csv(rmsf_path)
                
                pdb_data = {
                    "pdb_id": pdb,
                    "rmsf": rmsf_df
                }
                
                if os.path.exists(rmsf_bl_path):
                    bl_rmsf_df = pd.read_csv(rmsf_bl_path)
                    pdb_data["rmsf_bl"] = bl_rmsf_df
                
                pdb_rmsf_set.append(pdb_data)
            else:
                print(f"[rmsf_distribution.py] RMSF file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsf_distribution.py] Error processing {pdb}: {e}")

    make_dataset_rmsf_boxplot(pdb_rmsf_set, output_path, dataset_name, protein_class_mapping)