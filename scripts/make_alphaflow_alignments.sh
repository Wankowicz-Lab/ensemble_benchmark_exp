
if [ "$#" -ne 1 ]; then
    echo "[make_alphaflow_alignments.sh] A split name is required"
    exit 1
fi

SPLIT_NAME=$1
SPLIT_PATH="./splits/${SPLIT_NAME}_seqs.csv"

ABSOLUTE_SPLIT_PATH=$(realpath "$SPLIT_PATH")
ABSOLUTE_ALIGNMENTS_PATH=$(realpath "./bin/alignment/")
cd ./models/alphaflow/
python -m scripts.mmseqs_query --split $ABSOLUTE_SPLIT_PATH --outdir $ABSOLUTE_ALIGNMENTS_PATH
cd ../../
echo "[make_alphaflow_alignments.sh] Alignments for split ${SPLIT_NAME} are saved in ${ABSOLUTE_ALIGNMENTS_PATH}"