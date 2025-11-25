# python rmsf_cs_pdbs.py <dataset_name> <output_path>
# Makes a graph of the cosine similarity values with PDBs on x axis.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.lines import Line2D


def read_cosine_similarity(pdb_id):
    pdb_dir = f"./PDBs/{pdb_id}"
    cs_file = f"{pdb_dir}/analysis/cosine_similarity.csv"
    if os.path.exists(cs_file):
        return pd.read_csv(cs_file)
    else:
        print(f"[rmsf_cs_pdbs.py] Cosine similarity file not found for {pdb_id}")
        return None


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


def make_cs_scatter_plot(pdb_cs_list, output_path, dataset_name, protein_class_mapping=None):
    print(f"[rmsf_cs_pdbs.py] Making cosine similarity scatter plot for {len(pdb_cs_list)} PDBs...")
    
    if not pdb_cs_list:
        print("[rmsf_cs_pdbs.py] No data to plot")
        return

    # Define colors for the predictors
    predictorcolors = {
        'sam2': '#2ecc71',       # Green
        'alphaflow': '#e74c3c',   # Red
        'bioemu': '#3498db',      # Blue
        'openfold': '#f1c40f',    # Yellow
        'boltz2': '#9b59b6',      # Purple
    }

    predictors = ['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2']
    
    class_grouped_pdbs = {}
    if protein_class_mapping:
        pdb_to_class = {}
        for cls, pdbs in protein_class_mapping.items():
            for pdb in pdbs:
                pdb_to_class[pdb.lower()] = cls
        
        for pdb_data in pdb_cs_list:
            pdb_id = pdb_data["pdb_id"].lower()
            protein_class = pdb_to_class.get(pdb_id, "Other")
            
            if protein_class not in class_grouped_pdbs:
                class_grouped_pdbs[protein_class] = []
            class_grouped_pdbs[protein_class].append(pdb_data)
    else:
        class_grouped_pdbs["Other"] = pdb_cs_list
    
    fig, ax = plt.subplots(figsize=(24, 10))
    
    pdb_spacing = 3.0
    class_spacing = 0.0
    
    sorted_classes = ["Lysozyme", "CypA", "UBI", "STN", "Other"]
    
    current_pos = 1
    pdb_positions = {}
    class_boundaries = []
    class_centers = []
    class_names = []
    pdb_centers = []
    pdb_labels = []
    
    for class_idx, protein_class in enumerate(sorted_classes):
        if protein_class not in class_grouped_pdbs:
            continue
            
        pdbs_in_class = class_grouped_pdbs[protein_class]
        class_start_pos = current_pos
        
        for pdb_data in pdbs_in_class:
            pdb_id = pdb_data["pdb_id"]
            pdb_positions[pdb_id] = current_pos
            pdb_centers.append(current_pos)
            pdb_labels.append(pdb_id.upper())
            current_pos += pdb_spacing
        
        class_centers.append((class_start_pos + current_pos - pdb_spacing) / 2)
        class_names.append(protein_class)
        
        if class_idx < len(sorted_classes) - 1:
            class_boundaries.append(current_pos)
            current_pos += class_spacing
    
    marker_styles = ['o','o','o','o','o']#['o', 's', 'D', '^', 'v']  # Different marker styles
    marker_size = 150
    
    for class_idx, protein_class in enumerate(sorted_classes):
        if protein_class not in class_grouped_pdbs:
            continue
            
        pdbs_in_class = class_grouped_pdbs[protein_class]
        
        for pdb_data in pdbs_in_class:
            pdb_id = pdb_data["pdb_id"]
            cs_df = pdb_data["cosine_similarity"]
            
            pos = pdb_positions[pdb_id]
            
            for i, predictor in enumerate(predictors):
                cs_value = None
                
                predictor_rows = cs_df[
                    ((cs_df['element1'] == predictor) & (cs_df['element2'] == 'deposited')) |
                    ((cs_df['element1'] == 'deposited') & (cs_df['element2'] == predictor))
                ]
                
                if not predictor_rows.empty:
                    cs_value = predictor_rows['cosine_similarity'].values[0]
                    ax.scatter(pos, cs_value, color=predictorcolors[predictor], 
                               s=marker_size, marker=marker_styles[i], 
                               edgecolor='black', linewidth=1, alpha=0.8)
    
    for boundary in class_boundaries:
        ax.axvline(x=boundary - pdb_spacing/2, color='black', linestyle='--', alpha=0.5, linewidth=2)
    
    for center, class_name in zip(class_centers, class_names):
        ax.annotate(class_name, xy=(center, -0.08), xycoords=('data', 'axes fraction'),
                   ha='center', va='center', fontsize=24, fontweight='bold')
    
    ax.set_xticks(pdb_centers)
    ax.set_xticklabels(pdb_labels, fontsize=22, rotation=0)
    
    ax.set_ylim(0, 1.0)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.set_yticklabels([f"{y:.1f}" for y in np.arange(0, 1.1, 0.1)], fontsize=20)
    
    legend_elements = []
    for i, predictor in enumerate(predictors):
        legend_elements.append(
            Line2D([0], [0], marker=marker_styles[i], color='w', markerfacecolor=predictorcolors[predictor],
                   markersize=15, label=predictor.upper(), markeredgecolor='black')
        )
    
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.05), 
              fontsize=24, framealpha=0.9)
    
    ax.set_ylabel('Cosine Similarity', fontsize=36)
    #ax.set_title(f'Predictor-Deposited Cosine Similarity by PDB - {dataset_name}', fontsize=28)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    print(f"[rmsf_cs_pdbs.py] Cosine similarity scatter plot saved to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rmsf_cs_pdbs.py <dataset_name> <output_path>")
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
    
    pdb_cs_set = []
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
            cs_path = f"./PDBs/{pdb}/analysis/cosine_similarity.csv"
            
            if os.path.exists(cs_path):
                cs_df = pd.read_csv(cs_path)
                
                pdb_data = {
                    "pdb_id": pdb,
                    "cosine_similarity": cs_df
                }
                
                pdb_cs_set.append(pdb_data)
            else:
                print(f"[rmsf_cs_pdbs.py] Cosine similarity file not found for {pdb}, skipping...")
                
        except Exception as e:
            print(f"[rmsf_cs_pdbs.py] Error processing {pdb}: {e}")

    make_cs_scatter_plot(pdb_cs_set, output_path, dataset_name, protein_class_mapping)