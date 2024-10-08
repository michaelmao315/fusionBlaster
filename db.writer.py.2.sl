#!/bin/bash
#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36000

input=$1
genome=$2
OUTPUT_DIR=$3  # Pass the output directory as a third argument

# Check if output directory is provided; otherwise, set a default
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="database"  # Default directory name
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check genome and set transcriptome accordingly
if [[ $genome == "grch37" || $genome == "hg19" ]]; then
    transcriptome=/project/hlilab/sam/transcriptome/grch37.transcriptome.tsv
elif [[ $genome == "grch38" || $genome == "hg38" ]]; then
    transcriptome=/project/hlilab/sam/transcriptome/grch38.transcriptome.tsv
elif [[ $genome == "grcm38" || $genome == "mm10" ]]; then
    transcriptome=/project/hlilab/sam/transcriptome/grcm38.transcriptome.tsv
else
    echo -e "Please Try again with one of the supported genomes below:\n grch37 or hg19\n grch38 or hg38\n grcm38 or mm10\n"
    exit 1
fi

# Load the necessary modules
module load gcc
module load anaconda
module load goolf R

# Define output files in the output directory
full_output="$OUTPUT_DIR/full_output.tsv"
filtered_output="$OUTPUT_DIR/filtered_output.tsv"
filtered_input="$OUTPUT_DIR/filtered_${input##*/}"

# Run the R script and save the full output in the output directory
paste <(cat "$input") <(Rscript /project/hlilab/software/fusionBlaster/genome2transcriptome.R $genome <(cut -f2- "$input")) > "$full_output"

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
python /project/hlilab/software/fusionBlaster/db.writer.py "$OUTPUT_DIR/python_input.tsv" $genome > "$OUTPUT_DIR/refDB.fa"

# Clean up temporary files
# rm -f "$full_output" "python_input.tsv"