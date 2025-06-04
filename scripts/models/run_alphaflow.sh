#!/bin/bash
# run_alphaflow.sh <split_name> <sequence_count>

if [ "$#" -ne 2 ]; then
    echo "[run_alphaflow.sh] A pdb id AND sequence count is required."
    exit 1
fi

SPLIT_NAME=$1

SEQ_FILE="./splits/${SPLIT_NAME}_seqs.csv"
MSA_DIR="./bin/alignment"
WEIGHTS_PATH="./bin/weights/alphaflow_pdb_base_202402.pt"
ALPHAFLOW_OUT_DIR="./bin/alphaflow_out/${SPLIT_NAME}"

if [ ! -f "$SEQ_FILE" ]; then
    echo "Sequence file not found: $SEQ_FILE"
    exit 1
fi

ABSOLUTE_OUTPUT_DIR=$(realpath $ALPHAFLOW_OUT_DIR)
ABSOLUTE_MSA_DIR=$(realpath $MSA_DIR)
ABSOLUTE_WEIGHTS_PATH=$(realpath $WEIGHTS_PATH)

mkdir -p $ALPHAFLOW_OUT_DIR

cd ./models/alphaflow

python predict.py --mode alphafold --input_csv $SEQ_FILE --msa_dir $ABSOLUTE_MSA_DIR --weights $ABSOLUTE_WEIGHTS_PATH --samples $2 --outpdb $ABSOLUTE_OUTPUT_DIR

cd ../..