# python ediam_scatter.py <pdb_id> <predictor> <output_path>
# Creates a scatter plot of EDIA values by predictor for a single PDB

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import json
import seaborn as sns

predictors = ['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2']

def make_ediam_scatter(pdb_id, predictor, df, output_path):
    print(f"[ediam_scatter.py] Creating scatter plot for {pdb_id}...")
    
    if df.empty:
        print("[ediam_scatter.py] No data to plot")
        return
    df = df[df['predictor'] == predictor]
    
    if 'residue' in df.columns:
        df = df.sort_values('residue')
    
    plt.figure(figsize=(12, 6), dpi=100)
    
    ax = plt.scatter(df['residue'], df['EDIAm'], c=df['EDIAm'], cmap='coolwarm_r', 
                    alpha=0.8, s=30, edgecolors='none')
    
    
    plt.xlabel('Residue Number', fontsize=42)
    plt.ylabel('EDIAm Value', fontsize=42)
    plt.title(f'EDIAm Values by Residue for {pdb_id.upper()} ({predictor})', fontsize=28)

    plt.xticks(fontsize=24)
    
    plt.tight_layout()
    
    plt.savefig(output_path)
    print(f"[ediam_scatter.py] Scatter plot saved to {output_path}")
    
    pymol_commands = generate_pymol_commands(pdb_id, predictor, df)
    print(pymol_commands)
    

    
    
def generate_pymol_commands(pdb_id, predictor, df):
    #avg_df = df.groupby(['residue', 'chain', 'aa']).agg({'EDIAm': 'mean'}).reset_index()
    commands = f"""# PyMOL commands for visualizing EDIAm values for {pdb_id} ({predictor})
load PDBs/{pdb_id}/models/{predictor}.pdb, {predictor}

alter all, b=0.0

"""
    
    commands += "# Set b-factors to EDIAm values\n"
    for _, row in df.iterrows():
        if pd.notna(row['residue']) and pd.notna(row['EDIAm']) and pd.notna(row['chain']) and pd.notna(row['frame']):
            model_name = f"{predictor}_{1+row['frame']:04d}"
            commands += f"alter {model_name} and resi {int(row['residue'])}, b={row['EDIAm']}\n"
    
    commands += """
spectrum b, red_white_blue, minimum=0, maximum=1
show cartoon
"""
    
    with open('pymol_commands.txt', 'w') as f:
        f.write(commands)
    return commands


def parse_density_fitness_csv(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.json"
    
    if not os.path.exists(file_path):
        print(f"[ediam_scatter.py] File not found: {file_path}")
        return pd.DataFrame()
    
    try:
   
        all_data = []
        parsed = json.load(open(file_path, 'r'))
        predictors = parsed.get('predictor', [])
        frames = parsed.get('frame', [])
        metrics = parsed.get('metrics', [])
        if not predictors or not frames or not metrics:
            print(f"[ediam_scatter.py] No valid data found in {file_path}")
            return pd.DataFrame()
        
        allKeys = predictors.keys()
        for key in allKeys:
            predictor = predictors[key]
            frame = frames[key]
            thismetrics = metrics[key]
            for metric in thismetrics:
                all_data.append({
                    'predictor': predictor,
                    'frame': frame,
                    'residue': metric.get('seqID', None),
                    'aa': metric.get('compID', None),
                    'RSCCS': metric.get('RSCCS', None),
                    'RSR': metric.get('RSR', None),
                    'EDIAm': metric.get('EDIAm', None),
                    'chain': metric.get('asymID', None)
                })
      
        return pd.DataFrame(all_data)
            
    except Exception as e:
        print(f"[ediam_scatter.py] Error reading {file_path}: {e}")
        return pd.DataFrame()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python ediam_scatter.py <pdb_id> <predictor> <output_path>")
        sys.exit(1)

    pdb_id = sys.argv[1].lower()
    predictor = sys.argv[2]
    output_path = sys.argv[3]
    
    ediam_data = parse_density_fitness_csv(pdb_id)
    
    if ediam_data.empty:
        print(f"[ediam_scatter.py] No ediam data found for {pdb_id}")
        sys.exit(1)
    
    make_ediam_scatter(pdb_id, predictor, ediam_data, output_path)