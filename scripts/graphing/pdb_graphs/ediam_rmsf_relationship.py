import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import json
import seaborn as sns
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch


def make_rmsf_ediam_plot(pdb_id, rmsf_df, ediam_df, output_path, predictor):
    print(f"[ediam_rmsf_relationship.py] Creating plot for {pdb_id}...")
    
    if rmsf_df.empty or ediam_df.empty:
        print("[ediam_rmsf_relationship.py] Not enough data to plot")
        return
    
    plt.figure(figsize=(16, 10))
    sns.set_style("whitegrid")
    
    merged_data = pd.merge(
        rmsf_df,
        ediam_df,
        on=['predictor', 'residue'],
        how='inner',
        suffixes=('_rmsf', '_ediam')
    )
    
    if merged_data.empty:
        print("[ediam_rmsf_relationship.py] No matching residues between datasets")
        return
    
    pred_data = merged_data[merged_data['predictor'] == predictor]
    
    if not pred_data.empty:
        pred_data = pred_data.sort_values('residue')
        
        ax = plt.gca()
        
       
        plt.plot(
            pred_data['residue'],
            pred_data['rmsf'],
            linestyle='-',
            color='red',
            alpha=0.9,
            lw=2,
            label='RMSF (Flexibility)'
        )
        
        plt.plot(
            pred_data['residue'],
            -pred_data['EDIAm'],
            linestyle='-',
            color='purple',
            alpha=0.9,
            lw=2,
            label='EDIAm (Density Fit - inverted)'
        )
        
        plt.grid(True, linestyle='--', alpha=0.3, zorder=1)
        
        # Create combined legend
        data_legend_elements = [
            plt.Line2D([0], [0], color='red', lw=2, label='RMSF (Flexibility)'),
            plt.Line2D([0], [0], color='purple', lw=2, label='EDIAm (Density Fit - inverted)')
        ]

        plt.title(f"RMSF vs EDIAm for {pdb_id} - Predictor: {predictor}", fontsize=32)
        plt.legend(handles=data_legend_elements, loc='upper right', fontsize=24)
        
        plt.xlabel('Residue Number', fontsize=42)
        plt.ylabel('RMSF (Å) / -EDIAm', fontsize=42)

        
        ax.tick_params(axis='both', which='major', labelsize=24)

        plt.axhline(y=0, color='black', linestyle='-', alpha=0.5, lw=1.5, zorder=2)
        
        plt.text(0.02, 0.98, 'High Flexibility', transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='red', alpha=0.1))
        plt.text(0.02, 0.02, 'High Density Fit', transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='blue', alpha=0.1))

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"[ediam_rmsf_relationship.py] Saved plot to {output_path}")
        
    else:
        print(f"[ediam_rmsf_relationship.py] No data found for predictor {predictor}")


def load_rmsf_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/rmsf.csv"
    
    if not os.path.exists(file_path):
        print(f"[ediam_rmsf_relationship.py] RMSF file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"[ediam_rmsf_relationship.py] Error reading RMSF data: {e}")
        return pd.DataFrame()


def parse_density_fitness_csv(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.json"
    
    if not os.path.exists(file_path):
        print(f"[ediam_rmsf_relationship.py] File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        all_data = []
        parsed = json.load(open(file_path, 'r'))
        predictors = parsed.get('predictor', [])
        frames = parsed.get('frame', [])
        metrics = parsed.get('metrics', [])
        if not predictors or not frames or not metrics:
            print(f"[ediam_rmsf_relationship.py] No valid data found in {file_path}")
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
           
        ediam_df = pd.DataFrame(all_data)
            
        if not ediam_df.empty:
            avg_ediam = ediam_df.groupby(['predictor', 'residue', 'aa'])['EDIAm'].mean().reset_index()
            return avg_ediam
    except Exception as e:
        print(f"[ediam_rmsf_relationship.py] Error reading {file_path}: {e}")
        return pd.DataFrame()


def load_ediam_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.csv"
    
    if not os.path.exists(file_path):
        print(f"[ediam_rmsf_relationship.py] Density fitness file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        lines = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    lines.append({
                        'predictor': line.split(',')[0].strip(),
                        'frame': line.split(',')[1].strip(),
                        'metrics': (','.join(line.split(',')[2:]).strip())[1:-1]
                    })
        
        df = pd.DataFrame(lines)
        all_data = []
        
        if 'metrics' in df.columns:
            for _, row in df.iterrows():
                predictor = row['predictor']
                frame = row['frame']
                
                try:
                    metrics = json.loads(row['metrics'])
                    
                    for residue_data in metrics:
                        if 'EDIAm' in residue_data and 'seqID' in residue_data:
                            all_data.append({
                                'predictor': predictor,
                                'frame': frame,
                                'residue': residue_data['pdb']['seqNum'],
                                'EDIAm': residue_data['EDIAm'],
                                'aa': residue_data.get('compID', '')
                            })
                except Exception as e:
                    print(f"[ediam_rmsf_relationship.py] Error parsing metrics: {e}")
            
            ediam_df = pd.DataFrame(all_data)
            
            if not ediam_df.empty:
                avg_ediam = ediam_df.groupby(['predictor', 'residue', 'aa'])['EDIAm'].mean().reset_index()
                return avg_ediam
        
        return pd.DataFrame()
        
    except Exception as e:
        print(f"[ediam_rmsf_relationship.py] Error loading density fitness data: {e}")
        return pd.DataFrame()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python ediam_rmsf_relationship.py <pdb_id> <predictor> <output_path>")
        sys.exit(1)

    pdb_id = sys.argv[1].lower()
    predictor = sys.argv[2].lower()
    output_path = sys.argv[3]
    
    rmsf_data = load_rmsf_data(pdb_id)
    ediam_data = parse_density_fitness_csv(pdb_id)
    
    if rmsf_data.empty:
        print(f"[ediam_rmsf_relationship.py] No RMSF data found for {pdb_id} with predictor {predictor}")
        sys.exit(1)
        
    if ediam_data.empty:
        print(f"[ediam_rmsf_relationship.py] No EDIAm data found for {pdb_id} with predictor {predictor}")
        sys.exit(1)
    
    make_rmsf_ediam_plot(pdb_id, rmsf_data, ediam_data, output_path, predictor)