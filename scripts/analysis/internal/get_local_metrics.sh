#!/bin/bash
# get pdb local metrics using CCP4 density-fitness

MTZ_PATH=$1
PDB_PATH=$2
OUTPUT_PATH=$3

echo $PDB_PATH
echo $MTZ_PATH
echo $OUTPUT_PATH

if [ -z "$MTZ_PATH" ] || [ -z "$PDB_PATH" ] || [ -z "$OUTPUT_PATH" ]; then
    echo "Usage: $0 <mtz_path> <pdb_path> <output_path>"
    exit 1
fi

source /home/sbgrid/programs/x86_64-linux/ccp4/9.0/ccp4-9.0/bin/ccp4.setup-sh 
density-fitness "$MTZ_PATH" "$PDB_PATH" -o "$OUTPUT_PATH"