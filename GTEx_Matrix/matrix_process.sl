#!/bin/bash
#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=0-23:59:59
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=36000

# Check if all required arguments are provided
if [ $# -ne 5 ]; then
    echo "Usage: sbatch $0 <input_matrix> <sample_info_file> <output_dir> <group1_name> <group2_name>"
    exit 1
fi

INPUT_MATRIX="$1"
SAMPLE_INFO_FILE="$2"
OUTPUT_DIR="$3"
GROUP1_NAME="$4"
GROUP2_NAME="$5"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set predefined output file names
TTEST_OUTPUT="$OUTPUT_DIR/matrix_ttest.tsv"
CHI_SQUARE_OUTPUT="$OUTPUT_DIR/matrix_t_chi.tsv"
FINAL_OUTPUT="$OUTPUT_DIR/final.tsv"
SHORT_OUTPUT="$OUTPUT_DIR/short.tsv"

# Set fixed paths for the Python scripts
TTEST_CALCULATOR_SCRIPT="./fusionBlaster/GTEx_Matrix/calculate_T_test.py"
CHI_SQUARE_CALCULATOR_SCRIPT="./fusionBlaster/GTEx_Matrix/calculate_chi.py"
PROCESS_CHI_SQUARE_SCRIPT="./fusionBlaster/GTEx_Matrix/create_short.py"

# Load required modules
module load anaconda

# Run T-test calculator
echo "Running T-test calculation..."
python "$TTEST_CALCULATOR_SCRIPT" "$INPUT_MATRIX" "$SAMPLE_INFO_FILE" "$TTEST_OUTPUT" "$GROUP1_NAME" "$GROUP2_NAME"

# Run Chi-square calculator
echo "Running Chi-square calculation..."
python "$CHI_SQUARE_CALCULATOR_SCRIPT" "$TTEST_OUTPUT" "$SAMPLE_INFO_FILE" "$CHI_SQUARE_OUTPUT" "$GROUP1_NAME" "$GROUP2_NAME"

# Process Chi-square output
python "$PROCESS_CHI_SQUARE_SCRIPT" "$CHI_SQUARE_OUTPUT" "$FINAL_OUTPUT" "$SHORT_OUTPUT"

echo "Analysis complete. Final output is in $FINAL_OUTPUT"
echo "Short version is in $SHORT_OUTPUT"