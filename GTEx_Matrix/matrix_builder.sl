#!/bin/bash
#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36000

# Check if all arguments are provided
if [ $# -ne 5 ]; then
    echo "Usage: sbatch $0 <sample_sheet> <chimera_file> <refDB> <output_matrix> <level>"
    echo "  level: g for gene, b for breakpoint, t for transcript"
    exit 1
fi

sample_sheet=$1
chimera_file=$2
refDB=$3
output_matrix=$4
level=$5

# Set the appropriate script and input file based on the analysis level
case $level in
    g) 
        script="./fusionBlaster/GTEx_Matrix/gene_matrix_builder.py"
        input_file=$chimera_file
        ;;
    b) 
        script="./fusionBlaster/GTEx_Matrix/bp_matrix_builder.py"
        input_file=$chimera_file
        ;;
    t) 
        script="./fusionBlaster/GTEx_Matrix/transcript_matrix_builder.py"
        input_file=$refDB
        ;;
    *) 
        echo "Error: Invalid analysis level. Use g, b, or t."
        exit 1
        ;;
esac

# Load necessary modules (adjust as needed for your cluster)
module load anaconda

# Run the appropriate matrix builder script
python $script $sample_sheet $input_file $output_matrix

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Matrix building completed successfully. Output saved to $output_matrix"
else
    echo "Error: Matrix building failed. Check the error log for details."
    exit 1
fi