#!/bin/bash
# run_bioemu.sh <pdb_id> <sequence_count>

if [ "$#" -ne 2 ]; then
    echo "[run_bioemu.sh] A pdb id AND sequence count is required."
    exit 1
fi

PDB_ID=$1
PDB_DIR="./PDBs/${PDB_ID,,}"
SEQ_FILE="${PDB_DIR}/${PDB_ID,,}_seq.txt"
OUTPUT_FILE="${PDB_DIR}/${PDB_ID,,}_bioemu.pdb"

if [ ! -f "$SEQ_FILE" ]; then
    echo "Sequence file not found: $SEQ_FILE"
    exit 1
fi

SEQUENCE=$(cat "$SEQ_FILE")

export PYTHONPATH=$(pwd):$PYTHONPATH

python ./models/bioemu/inference.py "$PDB_ID" "$OUTPUT_FILE" "$SEQUENCE" $2
echo "[run_bioemu.sh] Saved $PDB_ID to $OUTPUT_FILE"