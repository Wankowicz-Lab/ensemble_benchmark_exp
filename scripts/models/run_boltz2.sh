#!/bin/bash
# run_boltz2.sh <pdb_id> <sequence_count>
set -e

if [ "$#" -ne 2 ]; then
    echo "[run_boltz2.sh] A pdb id AND sequence count is required."
    exit 1
fi

PDB_ID=$1
PDB_DIR="./PDBs/${PDB_ID,,}"
SEQ_FILE="${PDB_DIR}/${PDB_ID,,}_seq.txt"

if [ ! -f "$SEQ_FILE" ]; then
    echo "Sequence file not found: $SEQ_FILE"
    exit 1
fi


# Make fasta file needed by boltz2 before running
SEQUENCE=$(cat "$SEQ_FILE")
# YML_FILE="${PDB_DIR}/${PDB_ID,,}.fasta"

# if [ ! -f "$YML_FILE" ]; then
#     ALIGNMENT_LOCATION="./bin/alignment/$PDB_ID/a3m/$PDB_ID.a3m"
#     echo ">$PDB_ID|protein|$ALIGNMENT_LOCATION" > "$YML_FILE"
#     echo "$SEQUENCE" >> "$YML_FILE"
#     echo "[run_boltz2.sh] Created FASTA file for Boltz: $YML_FILE"
# fi


OUT_DIR="$PDB_DIR/boltz2_bin"
mkdir -p $OUT_DIR
rm -rf "$OUT_DIR/*"


YML_FILE="$OUT_DIR/boltz2.yml"
echo "version: 1" > "$YML_FILE"
echo "sequences:" >> "$YML_FILE"
echo "  - protein:" >> "$YML_FILE"
echo "     id: [A]" >> "$YML_FILE"
echo "     sequence: $SEQUENCE" >> "$YML_FILE"
echo "     msa: ./bin/alignment/$PDB_ID/a3m/$PDB_ID.a3m" >> "$YML_FILE"


# actually run boltz2

echo "[run_boltz2.sh] Running Boltz2 for $PDB_ID with sequence count $2"

boltz predict $YML_FILE --out_dir $OUT_DIR --output_format pdb --diffusion_samples $2 --override --use_msa_server

INPUT_FILE_BASENAME=$(basename "$YML_FILE" .yml)

PREDICTIONS_DIR="$OUT_DIR/boltz_results_$INPUT_FILE_BASENAME/predictions/$INPUT_FILE_BASENAME"
if [ ! -d "$PREDICTIONS_DIR" ]; then
    echo "[run_boltz2.sh] Unfortunately, No predictions found in $PREDICTIONS_DIR"
    exit 1
fi

# Compile all the prediction PDBs into a single PDB file with all conformations

echo "[run_boltz2.sh] Compiling predictions into a single PDB file"
OUTPUT_PDB="$PDB_DIR/boltz2_bin/${PDB_ID,,}_ensemble.pdb"

if [ -f "$OUTPUT_PDB" ]; then
    echo "[run_boltz2.sh] Removing existing output PDB file: $OUTPUT_PDB"
    rm "$OUTPUT_PDB"
fi

for PDB_FILE in "$PREDICTIONS_DIR"/*.pdb; do
    if [ -f "$PDB_FILE" ]; then
        echo "[run_boltz2.sh] Adding $PDB_FILE to $OUTPUT_PDB"
        echo "MODEL" >> "$OUTPUT_PDB"
        cat "$PDB_FILE" >> "$OUTPUT_PDB"
        echo "" >> "$OUTPUT_PDB"
        echo "ENDMDL" >> "$OUTPUT_PDB"
    fi
done

echo "[run_boltz2.sh] Boltz2 run completed for $PDB_ID. Aligning ensemble..."


# Align everything to the deposited PDB
python ./scripts/helpers/align_with_mdtraj.py "$OUTPUT_PDB" "$PDB_DIR/${PDB_ID,,}_final.pdb"


echo "[run_boltz2.sh] Alignment with MDTRAJ completed. Now aligning with PHASER..."
bash ./scripts/helpers/align_with_phaser.sh "${PDB_ID,,}" boltz2

echo "[run_boltz2.sh] Final, aligned PDB file: $PDB_DIR/${PDB_ID,,}_boltz2.pdb"