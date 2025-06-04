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

TOPOLOGY_FILE="${PDB_DIR}/${PDB_ID,,}_bioemu.pdb"
XTC_FILE="${PDB_DIR}/bioemu.xtc"


python -m bioemu.sidechain_relax --pdb-path $TOPOLOGY_FILE --xtc-path $XTC_FILE
echo "[run_bioemu.sh] Saved $PDB_ID to $OUTPUT_FILE"

#python ./scripts/helpers/align_to_exp_data.py "$OUTPUT_FILE" "$PDB_DIR/${PDB_ID,,}_final.pdb"

echo "[run_bioemu.sh] Aligned $OUTPUT_FILE to original"