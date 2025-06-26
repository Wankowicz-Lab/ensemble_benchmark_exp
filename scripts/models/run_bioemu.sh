#!/bin/bash
# run_bioemu.sh <pdb_id> <sequence_count>

set -e
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
rm -rf "$PDB_DIR/bioemu_bin/*"

XTC_FILE="${PDB_DIR}/bioemu_bin/bioemu.xtc"

python ./models/bioemu/inference.py "$PDB_ID" "$OUTPUT_FILE" "$SEQUENCE" $2
echo "[run_bioemu.sh] Saved backbone of $PDB_ID to $OUTPUT_FILE"

python -m bioemu.sidechain_relax --pdb-path $OUTPUT_FILE --xtc-path $XTC_FILE --outpath "$PDB_DIR/bioemu_bin"
echo "[run_bioemu.sh] Saved sidechains of $PDB_ID to $PDB_DIR/bioemu_bin"

ENSEMBLE_OUTPUT="${PDB_DIR}/bioemu_bin/${PDB_ID,,}_ensemble.pdb"

python ./scripts/helpers/dcd_to_pdb.py "$PDB_DIR/bioemu_bin/samples_sidechain_rec.xtc" "$PDB_DIR/bioemu_bin/samples_sidechain_rec.pdb" "$ENSEMBLE_OUTPUT"
echo "[run_bioemu.sh] Saved predicted ensemble of $PDB_ID to $ENSEMBLE_OUTPUT"

python ./scripts/helpers/align_with_mdtraj.py "$ENSEMBLE_OUTPUT" "$PDB_DIR/${PDB_ID,,}_final.pdb"
echo "[run_bioemu.sh] Aligned $ENSEMBLE_OUTPUT to original"