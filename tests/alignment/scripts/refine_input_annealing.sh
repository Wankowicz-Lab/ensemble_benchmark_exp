
set -e
#source /programs/sbgrid.shrc

PROTEIN_NAME="3k0n"
INPUT_DIR="./input/3k0n.openfold"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR not found"
    exit 1
fi


mkdir -p "$INPUT_DIR/phenix_refine_annealing"

MTZ_FILE=$(find "$INPUT_DIR" -name "${PROTEIN_NAME%.*}_final.mtz" 2>/dev/null | head -1)
    
if [ -z "$MTZ_FILE" ]; then
    echo "Warning: No MTZ file found for $PROTEIN_NAME, skipping..."
    continue
fi
    
echo "Using MTZ file: $MTZ_FILE"

for MODEL_PATH in "$INPUT_DIR/models"/*.pdb; do

    echo "========================================="
    echo "Processing: $MODEL_PATH"
    echo "========================================="

    MODEL_BASENAME=$(basename "$MODEL_PATH" .pdb)
    
    mkdir -p "$INPUT_DIR/models/phenix_refine_annealing_${MODEL_BASENAME}"


    ABS_MODEL_PATH=$(realpath "$MODEL_PATH")
    ABS_MTZ_PATH=$(realpath "$MTZ_FILE")

    cd "$INPUT_DIR/models/phenix_refine_annealing_${MODEL_BASENAME}"


    phenix.refine $ABS_MTZ_PATH $ABS_MODEL_PATH main.number_of_macro_cycles=8 strategy=rigid_body space_group="P 21 21 21" simulated_annealing=true
    cd - 
    cp "$INPUT_DIR/models/phenix_refine_annealing_${MODEL_BASENAME}/${MODEL_BASENAME}_refine_001.pdb" "$INPUT_DIR/phenix_refine_annealing/${MODEL_BASENAME}.pdb"
   

    echo "Finished processing: $MODEL_PATH"
done

OUTPUT_PDB="$INPUT_DIR/phenix_refined_annealing.pdb"
for REFINED_MODEL in "$INPUT_DIR/phenix_refine_annealing"/*.pdb; do
    if [ -f "$REFINED_MODEL" ]; then
        echo "Adding $REFINED_MODEL to $OUTPUT_PDB"
        echo "MODEL" >> "$OUTPUT_PDB"
        cat "$REFINED_MODEL" >> "$OUTPUT_PDB"
        echo "" >> "$OUTPUT_PDB"
        echo "ENDMDL" >> "$OUTPUT_PDB"
    fi
done
