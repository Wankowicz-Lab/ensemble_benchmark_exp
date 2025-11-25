
# if [ "$#" -ne 1 ]; then
#     echo "[make_alphaflow_alignments.sh] A split name is required"
#     exit 1
# fi

# SPLIT_NAME=$1
# SPLIT_PATH="./splits/${SPLIT_NAME}_seqs.csv"

# ABSOLUTE_SPLIT_PATH=$(realpath "$SPLIT_PATH")
# ABSOLUTE_ALIGNMENTS_PATH=$(realpath "./bin/alignment/")
# cd ./models/alphaflow/
# python -m scripts.mmseqs_query --split $ABSOLUTE_SPLIT_PATH --outdir $ABSOLUTE_ALIGNMENTS_PATH
# cd -
# echo "[make_alphaflow_alignments.sh] Alignments for split ${SPLIT_NAME} are saved in ${ABSOLUTE_ALIGNMENTS_PATH}"





### SPLIT INTO MULTI BATCHES VERSION ###

if [ "$#" -ne 1 ]; then
    echo "[make_alphaflow_alignments.sh] A split name is required"
    exit 1
fi

SPLIT_NAME=$1
SPLIT_PATH="./splits/${SPLIT_NAME}_seqs.csv"
BATCH_SIZE=15

ABSOLUTE_SPLIT_PATH=$(realpath "$SPLIT_PATH")
ABSOLUTE_ALIGNMENTS_PATH=$(realpath "./bin/alignment/")

TOTAL_LINES=$(tail -n +2 "$SPLIT_PATH" | wc -l)
NUM_BATCHES=$(( (TOTAL_LINES + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "[make_alphaflow_alignments.sh] Processing $TOTAL_LINES PDBs in $NUM_BATCHES batches of up to $BATCH_SIZE"

TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

head -n 1 "$SPLIT_PATH" > "$TEMP_DIR/header.csv"

for i in $(seq 0 $((NUM_BATCHES - 1))); do
    START_LINE=$((i * BATCH_SIZE + 2))
    END_LINE=$(((i + 1) * BATCH_SIZE + 1))
    
    BATCH_FILE="$TEMP_DIR/batch_${i}.csv"
    
    cat "$TEMP_DIR/header.csv" > "$BATCH_FILE"
    tail -n +2 "$SPLIT_PATH" | sed -n "${START_LINE},${END_LINE}p" >> "$BATCH_FILE"
    
    BATCH_COUNT=$(tail -n +2 "$BATCH_FILE" | wc -l)
    echo "[make_alphaflow_alignments.sh] Processing batch $((i + 1))/$NUM_BATCHES ($BATCH_COUNT PDBs)..."
    
    ABSOLUTE_BATCH_PATH=$(realpath "$BATCH_FILE")
    cd ./models/alphaflow/
    python -m scripts.mmseqs_query --split "$ABSOLUTE_BATCH_PATH" --outdir "$ABSOLUTE_ALIGNMENTS_PATH"
    cd -
done

echo "[make_alphaflow_alignments.sh] All alignments for split ${SPLIT_NAME} are saved in ${ABSOLUTE_ALIGNMENTS_PATH}"