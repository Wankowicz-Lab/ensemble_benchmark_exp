#!/bin/bash
# run_sam2.sh <pdb_id> <sequence_count>
set -e

if [ "$#" -ne 2 ]; then
    echo "[run_bioemu.sh] A pdb id AND sequence count is required."
    exit 1
fi

PDB_ID=$1
PDB_DIR="./PDBs/${PDB_ID,,}"
NOWAT_PDB_PATH="${PDB_DIR}/${PDB_ID,,}_nowat.pdb"
SEQ_FILE="${PDB_DIR}/${PDB_ID,,}_seq.txt"
OUTPUT_FILE="${PDB_DIR}/sam2_bin/${PDB_ID,,}_ensemble.pdb"

ABSOLUTE_PDB_DIR=$(realpath "$PDB_DIR")
ABSOLUTE_NOWAT_PDB_PATH=$(realpath "$NOWAT_PDB_PATH")
ABSOLUTE_SEQ_FILE=$(realpath "$SEQ_FILE")

# if [ -e "$OUTPUT_FILE" ]; then # can add this to avoid re-running if errors occur.
#     echo "[run_sam2.sh] Output file $OUTPUT_FILE already exists. Exiting."
#     exit 0
# fi

cd ./models/sam2 || exit 1

mkdir -p "$ABSOLUTE_PDB_DIR/sam2_bin"
rm -rf "$ABSOLUTE_PDB_DIR/sam2_bin/*"

python ./scripts/generate_ensemble.py -c ./config/mdcath_model.yaml -i $ABSOLUTE_NOWAT_PDB_PATH -o $ABSOLUTE_PDB_DIR/sam2_bin/sam2 -n $2 -b 8 -T 320 -d cuda
echo "[run_sam2.sh] Inference completed for $PDB_ID"

cd ../..

if [ ! -f "$ABSOLUTE_PDB_DIR/sam2_bin/sam2.traj.dcd" ]; then
    echo "[run_sam2.sh] Output file $ABSOLUTE_PDB_DIR/sam2_bin/sam2.traj.dcd does not exist. Exiting."
    exit 1
fi

python ./scripts/helpers/dcd_to_pdb.py "$ABSOLUTE_PDB_DIR/sam2_bin/sam2.traj.dcd" "$ABSOLUTE_PDB_DIR/sam2_bin/sam2.top.pdb" "$OUTPUT_FILE"
python ./scripts/helpers/align_with_mdtraj.py "$OUTPUT_FILE" "$PDB_DIR/${PDB_ID,,}_final.pdb"


echo "[run_sam2.sh] Alignment with MDTRAJ completed. Now aligning with PHASER..."
bash ./scripts/helpers/align_with_phaser.sh "${PDB_ID,,}" sam2
echo "[run_sam2.sh] Final, aligned PDB file: $PDB_DIR/${PDB_ID,,}_sam2.pdb"