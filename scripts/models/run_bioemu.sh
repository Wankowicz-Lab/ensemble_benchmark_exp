#!/bin/bash
# run_bioemu.sh <pdb_id> <sequence_count>

if [ "$#" -ne 2 ]; then
    echo "[run_bioemu.sh] A pdb id AND sequence count is required."
    exit 1
fi

PDB_ID=$1
PDB_DIR="./PDBs/${PDB_ID,,}"
SEQ_FILE="${PDB_DIR}/${PDB_ID,,}_seq.txt"
OUTPUT_FILE="${PDB_DIR}/bioemu_bin/backbone.pdb"

if [ ! -f "$SEQ_FILE" ]; then
    echo "Sequence file not found: $SEQ_FILE"
    exit 1
fi

SEQUENCE=$(cat "$SEQ_FILE")

export PYTHONPATH=$(pwd):$PYTHONPATH

mkdir -p "$PDB_DIR/bioemu_bin"

python ./models/bioemu/inference.py "$PDB_ID" "$OUTPUT_FILE" "$SEQUENCE" $2

XTC_FILE="${PDB_DIR}/bioemu.xtc"


python -m bioemu.sidechain_relax --pdb-path $OUTPUT_FILE --xtc-path $XTC_FILE --outpath "$PDB_DIR/bioemu_bin"
echo "[run_bioemu.sh] Saved backbone of $PDB_ID to $OUTPUT_FILE"

python ./scripts/helpers/align_to_exp_data.py "$PDB_DIR/bioemu_bin/samples_sidechain_rec.pdb" "$PDB_DIR/${PDB_ID,,}_final.pdb"
mv "$PDB_DIR/bioemu_bin/samples_sidechain_rec.pdb" "$PDB_DIR/${PDB_ID,,}_bioemu.pdb"

echo "[run_bioemu.sh] Aligned $OUTPUT_FILE to original"