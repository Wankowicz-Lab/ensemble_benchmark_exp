#!/bin/bash
# run_openfold.sh <split_name>
# note: openfold needs to be installed separately and the path will need to be changed:

set -e

OPENFOLD_PATH="/home/aslamaj/workdir/openfold/bin" # (contains openfold clone)


# Validation 

if [ ! -d "$OPENFOLD_PATH" ]; then
    echo "[run_openfold.sh] OpenFold path does not exist: $OPENFOLD_PATH. Modify the OPENFOLD_PATH variable in this script to the containing directory of openfold."
    exit 1
fi

if [ "$#" -ne 1 ]; then
    echo "[run_openfold.sh] A split name is required."
    exit 1
fi

SPLIT_NAME=$1
BIN_ROOT="./bin/openfold_data/$SPLIT_NAME"
ABS_BIN_ROOT=$(realpath $BIN_ROOT)

if [ ! -d "$ABS_BIN_ROOT" ]; then
    echo "[run_openfold.sh] The bin directory for the split does not exist: $ABS_BIN_ROOT. Please run the data preparation script first."
    exit 1
fi


# Get ready for output 
mkdir -p "$ABS_BIN_ROOT/outputs"
rm -rf "$ABS_BIN_ROOT/outputs/*"


echo "[run_openfold.sh] Running Openfold for $SPLIT_NAME"

cd $OPENFOLD_PATH || exit 1

RELRESOURCES="./openfold/resources"
export BASE_DATA_DIR=$(realpath $RELRESOURCES)

python ./openfold/run_pretrained_openfold.py \
    $ABS_BIN_ROOT/inputs \
    $ABS_BIN_ROOT/cifs \
    --output_dir $ABS_BIN_ROOT/outputs \
    --config_preset "model_3" \
    --model_device "cuda:0" \
    --use_precomputed_alignments $ABS_BIN_ROOT/alignments

cd - 
echo "[run_openfold.sh] OpenFold inference exited for $SPLIT_NAME. Processing outputs..."


# openfold.py never exits, so user has to manually terminate after completion. 

# Process outputs to our ./PDBs/*, and align each.

OUT_DIR="$ABS_BIN_ROOT/outputs/predictions"

for file in "$OUT_DIR"/*.pdb; do
    if [ -f "$file" ]; then
        if [[ ! "$file" =~ _model_3_relaxed\.pdb$ ]]; then
            echo "[run_openfold.sh] Skipping file: $file (not a model_3_relaxed output)"
            continue
        fi
        echo "[run_openfold.sh] Processing output file: $file"
        # 1pw2_model_3_relaxed.pdb -> 1pw2.pdb
        BASENAME=$(basename "$file")
        BASENAME_NO_EXT="${BASENAME%.*}"
        PDB_ID="${BASENAME_NO_EXT/_model_3_relaxed/}"
        
        echo $PDB_ID

        mkdir -p "./PDBs/$PDB_ID/openfold_bin"
        OUTPUT_PDB="./PDBs/$PDB_ID/openfold_bin/${PDB_ID}_ensemble.pdb" # NOT an ensemble, but consistency

        cp "$file" "$OUTPUT_PDB"
        echo "[run_openfold.sh] Copied $PDB_ID openfold output to $OUTPUT_PDB"

        python ./scripts/helpers/align_with_mdtraj.py "$OUTPUT_PDB" "./PDBs/${PDB_ID,,}/${PDB_ID,,}_final.pdb"

        echo "[run_openfold.sh] Alignment with MDTRAJ completed. Now aligning with PHASER..."
        bash ./scripts/helpers/align_with_phaser.sh "${PDB_ID,,}" openfold

        echo "[run_openfold.sh] Final, aligned PDB file: ./PDBs/${PDB_ID,,}/${PDB_ID,,}_openfold.pdb"

    fi
done

echo "[run_openfold.sh] OpenFold processing completed for $SPLIT_NAME."