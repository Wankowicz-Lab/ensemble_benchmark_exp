# python rscc_violin.py <pdb_id> <output_path>
# Creates a violin plot of RSCC values by predictor for a single PDB

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import json
import seaborn as sns

def make_rscc_violin(pdb_id, df, output_path):
    print(f"[rscc_violin.py] Creating violin plot for {pdb_id}...")
    
    if df.empty:
        print("[rscc_violin.py] No data to plot")
        return
    
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    
    colors = {
        'bioemu': '#3498db',    # Blue
        'alphaflow': '#e74c3c', # Red
        'sam2': '#2ecc71',       # Green
        'openfold': '#f1c40f',  # Yellow
        'boltz2': '#9b59b6',  # Purple
    }
    
    
    
    sns.stripplot(
        x='predictor', 
        y='RSCCS',
        data=df,
        palette=colors, 
        size=3,
        alpha=0.4, 
        jitter=True,
        dodge=False,
        order=['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2'], 
        zorder=1
    )
    
    ax = sns.violinplot(
        x='predictor', 
        y='RSCCS',
        data=df,
        palette=colors,
        inner='quartile',
        scale='width',
        order=['sam2', 'alphaflow', 'bioemu', 'openfold', 'boltz2'],
        alpha=0.75,
        zorder=10
    )
    
    means = df.groupby('predictor')['RSCCS'].mean()
    # for i, predictor in enumerate(['bioemu', 'alphaflow', 'sam2']):
    #     if predictor in means:
    #         plt.plot(i, means[predictor], 'o', color='white', markersize=10, markeredgecolor='black')
            
    plt.xlabel('Ensemble Predictor', fontsize=14)
    plt.ylabel('RSCCS (Real Space Correlation Coefficient)', fontsize=14)
    plt.title(f'RSCC of Residues of Conformations by Ensemble Predictions - {pdb_id.upper()}', fontsize=16)
    
    plt.ylim(-0.05, 1.05)
    
    plt.axhline(y=0.8, color="#376632", linestyle='--', alpha=0.5)
    #plt.text(plt.xlim()[0]-0.1, 0.8, '0.8 Good', va='center', ha='right', alpha=0.7)
    
    counts = df.groupby('predictor')['RSCCS'].count()
    ax.set_xticklabels([
        f"{predictor}\n(n={counts.get(predictor, 0)})"
        for predictor in ['sam2', 'alphaflow', 'bioemu', 
                          'openfold', 'boltz2']
    ])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"[rscc_violin.py] Saved violin plot to {output_path}")


def parse_density_fitness_csv(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.csv"
    
    if not os.path.exists(file_path):
        print(f"[rscc_violin.py] File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        lines = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    lines.append({
                        'predictor': line.split(',')[0].strip(),
                        'frame': line.split(',')[1].strip(),
                        'metrics': (','.join(line.split(',')[2:]).strip())[1:-1] # Remove surrounding quotes cuz JSON parsing
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
                        if 'RSCCS' in residue_data:
                            all_data.append({
                                'predictor': predictor,
                                'frame': frame,
                                'residue': residue_data.get('seqID', None),
                                'aa': residue_data.get('compID', None),
                                'RSCCS': residue_data['RSCCS'],
                                'RSR': residue_data.get('RSR', None),
                                'chain': residue_data.get('asymID', None)
                            })
                except json.JSONDecodeError:
                    print(row['metrics'])
                    print(f"[rscc_violin.py] Error parsing metrics JSON for {predictor} frame {frame}")
                except Exception as e:
                    print(f"[rscc_violin.py] Error processing metrics: {e}")
            
            return pd.DataFrame(all_data)
            
    except Exception as e:
        print(f"[rscc_violin.py] Error reading {file_path}: {e}")
        return pd.DataFrame()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rscc_violin.py <pdb_id> <output_path>")
        sys.exit(1)

    pdb_id = sys.argv[1].lower()
    output_path = sys.argv[2]
    
    rscc_data = parse_density_fitness_csv(pdb_id)
    
    if rscc_data.empty:
        print(f"[rscc_violin.py] No RSCC data found for {pdb_id}")
        sys.exit(1)
    
    make_rscc_violin(pdb_id, rscc_data, output_path)