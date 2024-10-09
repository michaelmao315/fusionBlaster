#!/bin/bash

input=$1
genome=$2
OUTPUT_DIR=$3  # Pass the output directory as a third argument

# Check if output directory is provided; otherwise, set a default
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="database"  # Default directory name
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check genome
case $genome in
    grch37|hg19|grch38|hg38|grcm38|mm10) ;;
    *)
        echo "Error: Unsupported genome '$genome'. Supported genomes are: grch37, hg19, grch38, hg38, grcm38, mm10" >&2
        exit 1
        ;;
esac

# Define output files in the output directory
full_output="$OUTPUT_DIR/full_output.tsv"
filtered_output="$OUTPUT_DIR/filtered_output.tsv"
filtered_input="$OUTPUT_DIR/filtered_${input##*/}"

# Run the R script and save the full output in the output directory
paste <(cat "$input") <(Rscript ./fusionBlaster/genome2transcriptome.R $genome <(cut -f2- "$input")) > "$full_output"

# Filter out intron regions and create filtered output and input files in the output directory
awk -F'\t' 'BEGIN {OFS="\t"} {
    if ($0 !~ /Intron Region/) {
        print $0 > "'"$filtered_output"'"
        print $1,$2,$3,$4,$5 > "'"$filtered_input"'"
    }
}' "$full_output"

# Extract only the R script output for the Python script
cut -f6- "$filtered_output" > "$OUTPUT_DIR/python_input.tsv"

# Call your Python script and save its output to refDB.fa in the output directory
python ./fusionBlaster/db.writer.py "$OUTPUT_DIR/python_input.tsv" $genome > "$OUTPUT_DIR/refDB.fa"

# Clean up temporary files
# rm -f "$full_output" "python_input.tsv"