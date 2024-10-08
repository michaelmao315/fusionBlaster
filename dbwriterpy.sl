#!/bin/bash
#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36000

input=$1
genome=$2

if [[ $genome == "grch37" || $genome == "hg19" ]]; then
        transcriptome=/project/hlilab/sam/transcriptome/grch37.transcriptome.tsv
elif [[ $genome == "grch38" || $genome == "hg38" ]]; then
        transcriptome=/project/hlilab/sam/transcriptome/grch38.transcriptome.tsv
elif [[ $genome == "grcm38" || $genome == "mm10" ]]; then
        transcriptome=/project/hlilab/sam/transcriptome/grcm38.transcriptome.tsv
else
        echo -e "Please Try again with one of the supported genomes below:\n grch37 or hg19\n grch38 or hg38\n grcm38 or mm10\n"
        exit
fi

# Load the necessary modules
module load gcc
module load anoconda 
module load goolf R

# Define a temporary file for intermediate output
temp_output=$(mktemp)

# Run the R script and filter the output, then write to a temporary file
Rscript ./fusionBlaster/genome2transcriptome.R $genome <(cut -f2- $input) | grep -v "Intron Region" > "$temp_output"

# Call your Python script with the temporary file as input
python ./fusionBlaster/db.writer.py "$temp_output" $genome

# Optionally, remove the temporary file if you don't need it after the script runs
rm "$temp_output"
