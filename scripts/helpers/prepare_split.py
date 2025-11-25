# prepare_split.py.py splits/dataset.txt
# Downloads and sets up ./PDBs/*

import os
import requests
import argparse
import pdb_sequence_maker
import gemmi

def get_rid_of_water_and_ions(pdb_file_path, output_file_path):
    cleaned_lines = []
    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if not line.startswith(('HETATM', 'ANISOU', 'TER', 'END')):
                cleaned_lines.append(line)
            # elif line.startswith('HETATM'):
            #     if not ('HOH' in line or 'WAT' in line or 'MG' in line):
            #         cleaned_lines.append(line)
    
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(cleaned_lines)
        output_file.close()
    return True

def download_file(url, output_path):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")
        return False


def download_pdb(pdb_id):
    pdb_id = pdb_id.lower()
    
    pdb_dir = os.path.join("PDBs", pdb_id)
    os.makedirs(pdb_dir, exist_ok=True)
    
    pdb_redo_url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb"
    pdb_output_path = os.path.join(pdb_dir, f"{pdb_id}_final.pdb")
    
    if os.path.exists(pdb_output_path):
        print(f"[prepare_split.py] PDB already exists: {pdb_output_path}")
        return True

    if download_file(pdb_redo_url, pdb_output_path):
        print(f"[prepare_split.py] Successfully downloaded PDB from PDB_REDO: {pdb_id}")
        
        mtz_url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.mtz"
        mtz_output_path = os.path.join(pdb_dir, f"{pdb_id}_final.mtz")
        
        if download_file(mtz_url, mtz_output_path):
            print(f"[prepare_split.py] Successfully downloaded MTZ from PDB_REDO: {pdb_id}")
        else:
            print(f"[prepare_split.py] Failed to download MTZ from PDB_REDO: {pdb_id}")

        cif_url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.cif"
        cif_output_path = os.path.join(pdb_dir, f"{pdb_id}_final.cif")
        
        if download_file(cif_url, cif_output_path):
            print(f"[prepare_split.py] Successfully downloaded CIF from PDB_REDO: {pdb_id}")
        else:
            print(f"[prepare_split.py] Failed to download CIF from PDB_REDO: {pdb_id}")

    else:
        print(f"[prepare_split.py] PDB_REDO not available for {pdb_id}. Trying RCSB...")
        
        rcsb_pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        rcsb_pdb_output = os.path.join(pdb_dir, f"{pdb_id}_final.pdb")
        
        if download_file(rcsb_pdb_url, rcsb_pdb_output):
            print(f"[prepare_split.py] Successfully downloaded PDB from RCSB: {pdb_id}")
            
            rcsb_cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            rcsb_cif_output = os.path.join(pdb_dir, f"{pdb_id}_final.cif")
            
            if download_file(rcsb_cif_url, rcsb_cif_output):
                print(f"[prepare_split.py] Successfully downloaded CIF from RCSB: {pdb_id}")
            else:
                print(f"[prepare_split.py] Failed to download CIF from RCSB: {pdb_id}")
            
            rcsb_sf_url = f"https://files.rcsb.org/download/{pdb_id}-sf.cif"
            rcsb_sf_output = os.path.join(pdb_dir, f"{pdb_id}_sf.cif")
            
            if download_file(rcsb_sf_url, rcsb_sf_output):
                print(f"[prepare_split.py] Successfully downloaded structure factors from RCSB: {pdb_id}")
                
                mtz_output = os.path.join(pdb_dir, f"{pdb_id}_final.mtz")
                if convert_cif_to_mtz(rcsb_sf_output, mtz_output, pdb_id):
                    os.rename(rcsb_sf_output, os.path.join(pdb_dir, f"{pdb_id}_sf_original.cif"))
                    
                else:
                    print(f"[prepare_split.py] Keeping structure factors in CIF format: {pdb_id}")
            else:
                print(f"[prepare_split.py] No structure factors available from RCSB: {pdb_id}")
                
        else:
            print(f"[prepare_split.py] Failed to download PDB for {pdb_id} from all sources")
            return False
    
    return True



def convert_cif_to_mtz(cif_path, mtz_path, pdb_id):
    try:
        sf_doc = gemmi.cif.read(cif_path)
        mtz = gemmi.Mtz()
        mtz.read_sf_mmcif(sf_doc)
        mtz.write_to_file(mtz_path)
        
        if os.path.exists(mtz_path):
            print(f"[prepare_split.py] Successfully converted CIF to MTZ using gemmi: {pdb_id}")
            return True
    except ImportError:
        print(f"[prepare_split.py] gemmi not available. Install with: pip install gemmi")
    except Exception as e:
        print(f"[prepare_split.py] gemmi conversion failed for {pdb_id}: {e}")
    
    print(f"[prepare_split.py] Could not convert CIF to MTZ for {pdb_id}.")
    return False


def main():
    parser = argparse.ArgumentParser(description='Initialize PDBs from a split file')
    parser.add_argument('split_file', help='Path to split text file')
    args = parser.parse_args()
    
  
    os.makedirs("PDBs", exist_ok=True)
    
    try:
        with open(args.split_file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"[prepare_split.py] Split file not found: {args.split_file}")
        return
    
    print(f"[prepare_split.py] Downloading {len(pdb_ids)} PDB structures...")
    
    current_iter = 0
    total_iter = len(pdb_ids)
    for pdb_id in pdb_ids:
        if download_pdb(pdb_id):
            # Downloaded the PDB successfully, now lets extract the sequence
            new_original_path = os.path.join("PDBs", pdb_id, f"{pdb_id}_final.pdb")
            thisPDBSeqArray = pdb_sequence_maker.get_sequence_array(new_original_path)
            if (len(thisPDBSeqArray) > 1):
                print(f"[prepare_split.py] Multiple sequences found for {pdb_id}, using the first one.")
            if len(thisPDBSeqArray) == 0:
                print(f"[prepare_split.py] No sequences found for {pdb_id}.")
                usingSequence = None
                continue
            usingSequence = thisPDBSeqArray[0]
            if usingSequence:
                seq_file_path = os.path.join("PDBs", pdb_id, f"{pdb_id}_seq.txt")
                with open(seq_file_path, 'w') as seq_file:
                    seq_file.write(usingSequence)
                print(f"[prepare_split.py] Sequence for {pdb_id} saved to {seq_file_path}")
            
            # Make non-water PDB
            cleaned_pdb_path = os.path.join("PDBs", pdb_id, f"{pdb_id}_nowat.pdb")
            if get_rid_of_water_and_ions(new_original_path, cleaned_pdb_path):
                print(f"[prepare_split.py] Cleaned PDB for {pdb_id} saved to {cleaned_pdb_path}")
            else:
                print(f"[prepare_split.py] Failed to clean PDB for {pdb_id}")
        current_iter += 1
        print(f"[prepare_split.py] Progress: {current_iter}/{total_iter} ({(current_iter / total_iter) * 100:.2f}%)")
    print("[prepare_split.py] All tasks completed.")

if __name__ == "__main__":
    main()