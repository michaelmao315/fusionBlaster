import sys
import pandas as pd

def combine_results(chimera_file, chi_square_file, output_file):
    # Read the chimera file
    chimera_headers = ['GeneSymbol', '5\'ENSG', '3\'ENSG', '5\'coord', '3\'coord']
    chimera_df = pd.read_csv(chimera_file, sep='\t', header=None, names=chimera_headers)
    
    print(f"Original chimera file lines: {len(chimera_df)}")
    
    # Create the connector column
    chimera_df['connector'] = chimera_df['5\'ENSG'] + 'x' + chimera_df['3\'ENSG']
    
    # Remove duplicates based on the connector column
    chimera_df = chimera_df.drop_duplicates(subset=['connector'], keep='first')
    
    print(f"Chimera file lines after removing duplicates: {len(chimera_df)}")
    
    # Read the Chi-square file
    chi_square_df = pd.read_csv(chi_square_file, sep='\t')
    
    print(f"Chi-square file lines: {len(chi_square_df)}")

    # Remove columns starting with 'SRR' from the Chi-square output
    chi_square_df = chi_square_df.loc[:, ~chi_square_df.columns.str.startswith('SRR')]

    # Merge the dataframes based on the connector column
    combined_df = pd.merge(chimera_df, chi_square_df, left_on='connector', right_on='Chimera', how='inner')

    # Remove the connector column as it's no longer needed
    combined_df = combined_df.drop('connector', axis=1)

    print(f"Combined dataframe lines: {len(combined_df)}")

    # Check if the number of lines in both original files is the same as the combined file
    if len(combined_df) != len(chi_square_df):
        print("Warning: The number of lines in the combined file does not match the original files.")
        print(f"Combined file: {len(combined_df)} lines")
        print(f"Original Chi-square file: {len(chi_square_df)} lines")

    # Save the combined result
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine_results.py <chimera_file> <chi_square_file> <output_file>")
        sys.exit(1)

    chimera_file = sys.argv[1]
    chi_square_file = sys.argv[2]
    output_file = sys.argv[3]

    combine_results(chimera_file, chi_square_file, output_file)