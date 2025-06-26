
if [ $# -eq 0 ]; then
  echo "Usage: $0 <path_to_split_file>"
  exit 1
fi

SPLIT_FILE=$1
if [ ! -f "$SPLIT_FILE" ]; then
  echo "Split file not found: $SPLIT_FILE"
  exit 1
fi

echo "[prepare_split.sh] Preparing split file: "$SPLIT_FILE""

echo "[prepare_split.sh] Normalizing split file: "$SPLIT_FILE""
contents=$(cat "$SPLIT_FILE")
normalized_contents=$(echo "$contents" | tr '[:upper:]' '[:lower:]')
echo "$normalized_contents" > "$SPLIT_FILE"

python ./scripts/helpers/prepare_split.py "$1"

SPLIT_BASENAME=$(basename "$SPLIT_FILE")
SPLIT_FILENAME_WITHOUT_EXT="${SPLIT_BASENAME%.*}"

echo "[prepare_split.sh] Generating split sequence CSV for: "$SPLIT_BASENAME""
python ./scripts/helpers/generate_split_seq_csv.py "$SPLIT_FILENAME_WITHOUT_EXT"

echo "[prepare_split.sh] Generating alphaflow alignments for: "$SPLIT_BASENAME""
bash "./scripts/helpers/make_alphaflow_alignments.sh" "$SPLIT_FILENAME_WITHOUT_EXT"

echo "[prepare_split.sh] Completed preparation for split file: "$SPLIT_BASENAME". It's ready to go with all models except openfold."

echo "[prepare_split.sh] Openfold's dir structure and alignments still need to be generated, which will search DBs and take some time. Use ./scripts/prepare_openfold.sh <split_name> to generate them."
