#!/bin/bash
# download_alignment_cifs.py <dataset>
# This script downloads the CIFs needed for OpenFold from the HHSearch hits

import os
import sys
import re
import urllib.request
import gzip
import shutil
import requests
import concurrent.futures
from threading import Lock


DATASET_ARG = sys.argv[1] if len(sys.argv) > 1 else None
if DATASET_ARG is None:
    print("Usage: python download_alignment_cifs.py <dataset>")
    sys.exit(1)

DATASET_ROOT = "./bin/openfold_data/" + DATASET_ARG;
ALIGNMENTS_DIR = DATASET_ROOT + "/alignments"
CIF_DIR = DATASET_ROOT + "/cifs"
INPUT_DIR = DATASET_ROOT + "/inputs"

print_lock = Lock()

def thread_safe_print(message):
    with print_lock:
        print(message)

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
        thread_safe_print(f"Warning: HHsearch file not found: {hhsearch_file}")
    return hits


def download_file(url, output_path):
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except requests.exceptions.RequestException as e:
        thread_safe_print(f"Error downloading {url}: {e}")
        return False


def download_cif(pdb_id, output_dir):
    cif_url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.cif"
    output_file = os.path.join(output_dir, f"{pdb_id}.cif")
    
    if os.path.exists(output_file):
        thread_safe_print(f"CIF already exists: {output_file}")
        return True
    
    try:
        thread_safe_print(f"Downloading {pdb_id}...")
        if download_file(cif_url, output_file):
            thread_safe_print(f"Downloaded {pdb_id} CIF to {output_file}")
            return True
        else:
            thread_safe_print(f"Failed to download {pdb_id} CIF")
            return False
    except Exception as e:
        thread_safe_print(f"Error downloading CIF for {pdb_id}: {e}")
        return False


def download_cifs_for_target(pdb_id, max_downloads=50, max_workers=10):
    alignment_hhsearch = f"{ALIGNMENTS_DIR}/{pdb_id.lower()}/hhsearch_output.hhr"
    
    hits = parse_hhsearch_hits(alignment_hhsearch)
    
    if not hits:
        thread_safe_print(f"No hits found for {pdb_id}")
        return
    
    thread_safe_print(f"Found {len(hits)} hits for {pdb_id}")
    
    hits_to_download = hits[:max_downloads]
    
    successful_downloads = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pdb = {
            executor.submit(download_cif, hit_pdb_id, CIF_DIR): hit_pdb_id 
            for hit_pdb_id in hits_to_download
        }
        
        for future in concurrent.futures.as_completed(future_to_pdb):
            hit_pdb_id = future_to_pdb[future]
            try:
                if future.result():
                    successful_downloads += 1
            except Exception as exc:
                thread_safe_print(f"Download generated an exception for {hit_pdb_id}: {exc}")
    
    thread_safe_print(f"Completed downloads for {pdb_id}: {successful_downloads}/{len(hits_to_download)} successful")


all_input_fastas = os.listdir(INPUT_DIR)
all_input_fastas = [f for f in all_input_fastas if f.endswith(".fasta")]
all_pdb_ids = [f.replace(".fasta", "") for f in all_input_fastas]

os.makedirs(CIF_DIR, exist_ok=True)

allCount = len(all_pdb_ids)
currCount = 0;
for pdb_id in all_pdb_ids:
    download_cifs_for_target(pdb_id, max_downloads=50, max_workers=10)
    currCount += 1
    print(f"Progress: {currCount}/{allCount} targets processed ({(currCount / allCount) * 100:.2f}%)")

print("CIF download process completed.")