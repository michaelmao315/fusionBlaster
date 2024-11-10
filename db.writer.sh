#!/bin/bash

# Get the directory where the script is located, resolving symlinks
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
FUSION_BLASTER_DIR="$(dirname "$SCRIPT_DIR")"

input=$1
genome=$2
OUTPUT_DIR=$3 # Pass the output directory as a third argument

# Check if output directory is provided; otherwise, set a default
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${SCRIPT_DIR}/database" # Default directory is now relative to script location
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check genome
case $genome in
    hg19)
        genome="grch37"
        ;;
    hg38)
        genome="grch38"
        ;;
    grch37|grch38) ;;
    *)
        echo "Error: Unsupported genome '$genome'. Supported genomes are: grch37, hg19, grch38, hg38" >&2
        exit 1
        ;;
esac

# Define output files in the output directory
full_output="$OUTPUT_DIR/full_output.tsv"
filtered_output="$OUTPUT_DIR/filtered_output.tsv"
filtered_input="$OUTPUT_DIR/filtered_${input##*/}"
python_input="$OUTPUT_DIR/python_input.tsv"

# Run the R script and save the full output in the output directory
paste <(cat "$input") <(Rscript "${SCRIPT_DIR}/genome2transcriptome.R" "$genome" <(cut -f2- "$input")) > "$full_output"

# Filter out intron regions and create filtered output and input files in the output directory
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($0 !~ /Intron Region/) {
        print $0 > "'"$filtered_output"'"
        print $1,$2,$3,$4,$5 > "'"$filtered_input"'"
    }
}' "$full_output"

# Extract only the R script output for the Python script
cut -f6- "$filtered_output" > "$python_input"

# Call your Python script and save its output to refDB.fa in the output directory
python "${SCRIPT_DIR}/db.writer.py" "$python_input" "$genome" "${SCRIPT_DIR}" > "$OUTPUT_DIR/refDB.fa"

# Clean up temporary files
# rm -f "$full_output" "$python_input"