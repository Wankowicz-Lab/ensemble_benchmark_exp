#!/bin/bash
# run_openfold.sh <split_name> [num_samples]
# note: openfold needs to be installed separately and the path will need to be changed:

#set -e

OPENFOLD_PATH="./models/openfold" # (contains openfold clone)



if [ ! -d "$OPENFOLD_PATH" ]; then
    echo "[run_openfold.sh] OpenFold path does not exist: $OPENFOLD_PATH. Modify the OPENFOLD_PATH variable in this script to the containing directory of openfold."
    exit 1
fi

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "[run_openfold.sh] Usage: $0 <split_name> [num_samples]"
    echo "[run_openfold.sh] Default num_samples is 5 if not specified."
    exit 1
fi


SPLIT_NAME=$1
NUM_SAMPLES=${2:-1} 

echo "[run_openfold.sh] Will generate $NUM_SAMPLES samples for each protein"

BIN_ROOT="./bin/openfold_data/$SPLIT_NAME"
ABS_BIN_ROOT=$(realpath $BIN_ROOT)

if [ ! -d "$ABS_BIN_ROOT" ]; then
    echo "[run_openfold.sh] The bin directory for the split does not exist: $ABS_BIN_ROOT. Please run the data preparation script first."
    exit 1
fi



########################
mkdir -p "$ABS_BIN_ROOT/outputs"
rm -rf "$ABS_BIN_ROOT/outputs/*"

echo "[run_openfold.sh] Running Openfold $NUM_SAMPLES times for $SPLIT_NAME"

cd $OPENFOLD_PATH || exit 1

RELRESOURCES="./openfold/resources"
export BASE_DATA_DIR=$(realpath $RELRESOURCES)

echo $(realpath $RELRESOURCES)

for i in $(seq 1 $NUM_SAMPLES); do
    echo "[run_openfold.sh] Running OpenFold sample $i/$NUM_SAMPLES..."
    
    SAMPLE_OUTPUT_DIR="$ABS_BIN_ROOT/outputs/sample_$i"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    # OpenFold hangs on completion so we're killing it after it finishes.

    python ./run_pretrained_openfold.py \
        $ABS_BIN_ROOT/inputs \
        $ABS_BIN_ROOT/cifs \
        --output_dir $SAMPLE_OUTPUT_DIR \
        --config_preset "model_3" \
        --model_device "cuda:0" \
        --use_precomputed_alignments $ABS_BIN_ROOT/alignments 2>&1 | \
    while IFS= read -r line; do
        echo "$line"
        
        if [[ "$line" =~ "please install the libaio" ]]; then # this string is outputted after openfold is done.
            echo "[run_openfold.sh] Detected completion signal, terminating OpenFold..."
            pkill -f "run_pretrained_openfold.py" || true
            break
        fi
        
        if [[ "$line" =~ "Error"|"Failed"|"Exception" ]]; then
            echo "[run_openfold.sh] Detected error, terminating OpenFold..."
            pkill -f "run_pretrained_openfold.py" || true
            break
        fi
    done
    
    sleep 1
    
    
    echo "[run_openfold.sh] Completed sample $i/$NUM_SAMPLES"
done

cd - 
#######################
echo "[run_openfold.sh] OpenFold inference completed for all $NUM_SAMPLES samples. Processing outputs..."

# Process outputs to our ./PDBs/*, concatenate samples, and align each.

FIRST_SAMPLE_DIR="$ABS_BIN_ROOT/outputs/sample_1/predictions"
for file in "$FIRST_SAMPLE_DIR"/*.pdb; do
    if [ -f "$file" ]; then
        if [[ ! "$file" =~ _model_3_relaxed\.pdb$ ]]; then
            echo "[run_openfold.sh] Skipping file: $file (not a model_3_relaxed output)"
            continue
        fi
        
        BASENAME=$(basename "$file")
        BASENAME_NO_EXT="${BASENAME%.*}"
        PDB_ID="${BASENAME_NO_EXT/_model_3_relaxed/}"
        
        echo "[run_openfold.sh] Processing protein: $PDB_ID"
        
        mkdir -p "./PDBs/$PDB_ID/openfold_bin"
        ENSEMBLE_PDB="./PDBs/$PDB_ID/openfold_bin/${PDB_ID}_ensemble.pdb"
        
        rm -f "$ENSEMBLE_PDB"
        
        MODEL_COUNT=1
        for i in $(seq 1 $NUM_SAMPLES); do
            SAMPLE_FILE="$ABS_BIN_ROOT/outputs/sample_$i/predictions/${PDB_ID}_model_3_relaxed.pdb"
            
            if [ -f "$SAMPLE_FILE" ]; then
                echo "[run_openfold.sh] Adding sample $i to ensemble for $PDB_ID"
                
                echo "MODEL     $MODEL_COUNT" >> "$ENSEMBLE_PDB"
                
                grep -v "^MODEL\|^ENDMDL" "$SAMPLE_FILE" >> "$ENSEMBLE_PDB"
                
                echo "ENDMDL" >> "$ENSEMBLE_PDB"
                
                MODEL_COUNT=$((MODEL_COUNT + 1))
            else
                echo "[run_openfold.sh] Warning: Sample file not found: $SAMPLE_FILE"
            fi
        done
        
        echo "[run_openfold.sh] Created ensemble PDB with $((MODEL_COUNT - 1)) models: $ENSEMBLE_PDB"
        
        # Align the ensemble
        python ./scripts/helpers/align_with_mdtraj.py "$ENSEMBLE_PDB" "./PDBs/${PDB_ID,,}/${PDB_ID,,}_final.pdb"

        echo "[run_openfold.sh] Alignment with MDTRAJ completed. Now aligning with PHASER..."
        bash ./scripts/helpers/align_with_phaser.sh "${PDB_ID,,}" openfold

        echo "[run_openfold.sh] Final, aligned PDB file: ./PDBs/${PDB_ID,,}/${PDB_ID,,}_openfold.pdb"
    fi
done

echo "[run_openfold.sh] OpenFold processing completed for $SPLIT_NAME with $NUM_SAMPLES samples per protein."