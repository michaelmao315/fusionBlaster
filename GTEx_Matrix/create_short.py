import sys
import pandas as pd

def process_chi_square_output(input_file, full_output_file, short_output_file):
    # Read the Chi-square file
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Original file lines: {len(df)}")
    
    # Remove columns starting with 'SRR'
    full_df = df.loc[:, ~df.columns.str.startswith('SRR')]
    
    # Save the full result (without SRR columns)
    full_df.to_csv(full_output_file, sep='\t', index=False)
    print(f"Full results (without SRR columns) saved to {full_output_file}")
    print(f"Number of lines in full output: {len(full_df)}")

    # Create short version
    columns_to_keep = ['Chimera', 'ChiSquare_Corrected_Pvalue', 'T_Corrected_Pvalue_5p', 'T_Corrected_Pvalue_3p']
    short_df = full_df[columns_to_keep]
    
    # Save the short version without header
    short_df.to_csv(short_output_file, sep='\t', index=False, header=False)
    print(f"Short version saved to {short_output_file}")
    print(f"Number of lines in short output: {len(short_df)}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_chi_square.py <input_file> <full_output_file> <short_output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    full_output_file = sys.argv[2]
    short_output_file = sys.argv[3]
    
    process_chi_square_output(input_file, full_output_file, short_output_file)