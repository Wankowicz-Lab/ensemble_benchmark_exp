#!/bin/bash
# run_alphaflow.sh <split_name> <sequence_count>
set -e

if [ "$#" -ne 2 ]; then
    echo "[run_alphaflow.sh] A split name AND sequence count is required."
    exit 1
fi


SPLIT_NAME=$1


SEQ_FILE="./splits/${SPLIT_NAME}_seqs.csv"
MSA_DIR="./bin/alignment"
WEIGHTS_PATH="./bin/weights/alphaflow_pdb_base_202402.pt"
ALPHAFLOW_OUT_DIR="./bin/alphaflow_out/${SPLIT_NAME}"

mkdir -p $ALPHAFLOW_OUT_DIR

if [ ! -f "$SEQ_FILE" ]; then
    echo "Sequence file not found: $SEQ_FILE"
    exit 1
fi


ABSOLUTE_OUTPUT_DIR=$(realpath $ALPHAFLOW_OUT_DIR)
ABSOLUTE_MSA_DIR=$(realpath $MSA_DIR)
ABSOLUTE_WEIGHTS_PATH=$(realpath $WEIGHTS_PATH)
ABSOLUTE_SEQ_FILE=$(realpath $SEQ_FILE)


mkdir -p $ALPHAFLOW_OUT_DIR
rm -rf $ALPHAFLOW_OUT_DIR
mkdir -p $ALPHAFLOW_OUT_DIR


echo "[run_alphaflow.sh] Starting AlphaFlow prediction for split: $SPLIT_NAME"



cd ./models/alphaflow


python predict.py --mode alphafold --input_csv $ABSOLUTE_SEQ_FILE --msa_dir $ABSOLUTE_MSA_DIR --weights $ABSOLUTE_WEIGHTS_PATH --samples $2 --outpdb $ABSOLUTE_OUTPUT_DIR

RUN_STATUS=$?
cd ../..

if [ $RUN_STATUS -eq 0 ]; then
    echo "[run_alphaflow.sh] AlphaFlow prediction completed successfully"
    
    echo "[run_alphaflow.sh] Moving PDB files from the output..."
    
    mkdir -p ./PDBs
    
    for pdb_file in $ALPHAFLOW_OUT_DIR/*.pdb; do
        if [ -f "$pdb_file" ]; then
            pdb_id=$(basename "$pdb_file" .pdb)
            
            target_dir="./PDBs/$pdb_id"
            mkdir -p "$target_dir"
            mkdir -p "$target_dir/alphaflow_bin"
            rm -rf "$target_dir/alphaflow_bin/*"
            OUTPUT_FILE="$target_dir/alphaflow_bin/${pdb_id}_ensemble.pdb"
            cp "$pdb_file" $OUTPUT_FILE
            python ./scripts/helpers/align_with_mdtraj.py $OUTPUT_FILE "$target_dir/${pdb_id,,}_final.pdb"

            echo "[run_alphaflow.sh] Alignment with MDTRAJ completed. Now aligning with PHASER..."
            bash ./scripts/helpers/align_with_phaser.sh "${pdb_id,,}" alphaflow
            echo "[run_alphaflow.sh] Final, aligned PDB file: $PDB_DIR/${pdb_id,,}_alphaflow.pdb"

        fi
    done
    
    echo "[run_alphaflow.sh] All AlphaFlow PDB files have been moved."
else
    echo "[run_alphaflow.sh] AlphaFlow prediction failed with status code $RUN_STATUS"
    exit $RUN_STATUS
fi

