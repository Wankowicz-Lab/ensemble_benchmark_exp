#!/usr/bin/env python3
# download_alignment_cifs.py <dataset>
# This script downloads the CIFs needed for OpenFold from the HHSearch hits using batch downloading

import os
import sys
import re
import subprocess
import tempfile
from pathlib import Path


DATASET_ARG = sys.argv[1] if len(sys.argv) > 1 else None
if DATASET_ARG is None:
    print("Usage: python download_alignment_cifs.py <dataset>")
    sys.exit(1)

DATASET_ROOT = "./bin/openfold_data/" + DATASET_ARG
ALIGNMENTS_DIR = DATASET_ROOT + "/alignments"
CIF_DIR = DATASET_ROOT + "/cifs"
INPUT_DIR = DATASET_ROOT + "/inputs"

def parse_hhsearch_hits(hhsearch_file):
    hits = []
    try:
        with open(hhsearch_file, 'r') as f:
            for line in f:
                if re.match(r'^\s*\d+\s+', line):
                    parts = line.split()
                    if len(parts) > 1:
                        pdb_entry = parts[1]
                        pdb_id = pdb_entry[:4].lower()
                        if len(pdb_id) == 4 and pdb_id.isalnum():
                            hits.append(pdb_id)
    except FileNotFoundError:
        print(f"Warning: HHsearch file not found: {hhsearch_file}")
    return hits


def collect_all_hits(all_pdb_ids, max_downloads=50):
    all_hits = set()
    
    for pdb_id in all_pdb_ids:
        alignment_hhsearch = f"{ALIGNMENTS_DIR}/{pdb_id.lower()}/hhsearch_output.hhr"
        hits = parse_hhsearch_hits(alignment_hhsearch)
        
        if not hits:
            print(f"No hits found for {pdb_id}")
            continue
        
        print(f"Found {len(hits)} hits for {pdb_id}")
        hits_to_add = hits[:max_downloads]
        all_hits.update(hits_to_add)
    
    return list(all_hits)


def batch_download_cifs(pdb_ids, output_dir):
    if not pdb_ids:
        print("No PDB IDs to download")
        return
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
        temp_file.write(','.join(pdb_ids))
        temp_file_path = temp_file.name
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        batch_script = "./scripts/helpers/pdbredo_batch_dl.sh"
        cmd = [
            "bash", batch_script,
            "-f", temp_file_path,
            "-o", output_dir,
            "-c" 
        ]
        
        print(f"Running batch download for {len(pdb_ids)} PDB IDs...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("Batch download completed successfully")
            print(result.stdout)
        else:
            print(f"Batch download failed with return code {result.returncode}")
            print(f"Error: {result.stderr}")
            
    finally:
        os.unlink(temp_file_path)


def rename_downloaded_files(output_dir):
    cif_dir = Path(output_dir)
    
    for gz_file in cif_dir.glob("*.cif.gz"):
        pdb_id = gz_file.stem.replace(".cif", "")
        target_file = cif_dir / f"{pdb_id}.cif"
        
        if target_file.exists():
            print(f"CIF already exists: {target_file}")
            gz_file.unlink()
            continue
        
        import gzip
        try:
            with gzip.open(gz_file, 'rt') as f_in:
                with open(target_file, 'w') as f_out:
                    f_out.write(f_in.read())
            
            print(f"Decompressed {gz_file} to {target_file}")
            gz_file.unlink()
            
        except Exception as e:
            print(f"Error decompressing {gz_file}: {e}")


def download_cifs_for_dataset():
    all_input_fastas = os.listdir(INPUT_DIR)
    all_input_fastas = [f for f in all_input_fastas if f.endswith(".fasta")]
    all_pdb_ids = [f.replace(".fasta", "") for f in all_input_fastas]
    
    os.makedirs(CIF_DIR, exist_ok=True)
    
    print(f"Processing {len(all_pdb_ids)} targets...")
    
    all_hits = collect_all_hits(all_pdb_ids, max_downloads=250)
    
    print(f"Found {len(all_hits)} unique PDB IDs to download")
    
    existing_cifs = set()
    for cif_file in Path(CIF_DIR).glob("*.cif"):
        pdb_id = cif_file.stem
        existing_cifs.add(pdb_id)
    
    hits_to_download = [hit for hit in all_hits if hit not in existing_cifs]
    
    if not hits_to_download:
        print("All required CIF files already exist!")
        return
    
    print(f"Need to download {len(hits_to_download)} new CIF files")
    
    batch_size = 100
    for i in range(0, len(hits_to_download), batch_size):
        batch = hits_to_download[i:i+batch_size]
        print(f"Downloading batch {i//batch_size + 1}/{(len(hits_to_download)-1)//batch_size + 1}")
        
        batch_download_cifs(batch, CIF_DIR)
        rename_downloaded_files(CIF_DIR)
    
    print("CIF download process completed.")


if __name__ == "__main__":
    download_cifs_for_dataset()
