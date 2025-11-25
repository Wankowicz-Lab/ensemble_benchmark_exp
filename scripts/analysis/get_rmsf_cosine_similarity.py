# python get_rmsf_cosine_similarity.py <pdb_id>
# adds a cosine_similarity.csv to PDB analysis folder

import sys
import os 
import pandas as pd
from itertools import combinations



# vector dimension is residues
# < RMSF, RMSF, RMSF, ... >
def get_cosine_similarity(vector1, vector2):
    dot_product = sum(a * b for a, b in zip(vector1, vector2))
    magnitude1 = sum(a * a for a in vector1) ** 0.5
    magnitude2 = sum(b * b for b in vector2) ** 0.5
    if magnitude1 == 0 or magnitude2 == 0:
        return 0.0
    return dot_product / (magnitude1 * magnitude2)



def get_predictor_rmsf(pdb):
    rmsf_file_path = f"./PDBs/{pdb}/analysis/rmsf.csv" # predictor,residue,residue_aa,rmsf
    df = pd.read_csv(rmsf_file_path)
    return df

def get_deposited_rmsf(pdb):
    rmsf_file_path = f"./PDBs/{pdb}/analysis/{pdb}_qfit_RMSF.csv" # ,resseq,AA,Chain,RMSF,PDB_name
    df = pd.read_csv(rmsf_file_path)
    return df


def get_rmsf_vectors(pdb):
    deposited_df = get_deposited_rmsf(pdb)
    predictor_df = get_predictor_rmsf(pdb)
    predictors = predictor_df['predictor'].unique()

    rmsf_vectors = {}
    rmsf_vectors['deposited'] = deposited_df['RMSF'].tolist()

    for predictor in predictors:
        predictor_data = predictor_df[predictor_df['predictor'] == predictor]
        rmsf_vectors[predictor] = predictor_data['rmsf'].tolist()
    
    all_elements = list(rmsf_vectors.keys()) # combination of ALL

    combination_data = []

    for elem1, elem2 in combinations(all_elements, 2):
        vec1 = rmsf_vectors[elem1]
        vec2 = rmsf_vectors[elem2]
        cos_sim = get_cosine_similarity(vec1, vec2)
        #print(f"Cosine similarity between {elem1} and {elem2}: {cos_sim:.4f}")
        combination_data.append((elem1, elem2, cos_sim))

    output_df = pd.DataFrame(combination_data, columns=['element1', 'element2', 'cosine_similarity'])
    output_path = f"./PDBs/{pdb}/analysis/cosine_similarity.csv"
    output_df.to_csv(output_path, index=False)
    print(f"[get_rmsf_cosine_similarity.py] Cosine similarity results saved to {output_path}")
    
    
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python get_rmsf_cosine_similarity.py <pdb_id>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    get_rmsf_vectors(pdb_id)