import sys
import pandas as pd
import numpy as np
from scipy import stats

def load_sample_info(sample_info_file):
    sample_info = {}
    with open(sample_info_file, 'r') as f:
        for line in f:
            srr, value, _ = line.strip().split('\t')
            sample_info[srr] = float(value)  # Convert to float for numerical analysis
    return sample_info

def calculate_correlations(matrix_df, sample_info):
    correlations_5p = []
    pvalues_5p = []
    correlations_3p = []
    pvalues_3p = []

    for _, row in matrix_df.iterrows():
        values = []
        ratios_5p = []
        ratios_3p = []

        for sample in sample_info.keys():
            if f"{sample}_5_ratio" in row and f"{sample}_3_ratio" in row:
                values.append(sample_info[sample])
                ratios_5p.append(float(row[f"{sample}_5_ratio"]))
                ratios_3p.append(float(row[f"{sample}_3_ratio"]))

        corr_5p, pval_5p = stats.pearsonr(values, ratios_5p)
        corr_3p, pval_3p = stats.pearsonr(values, ratios_3p)

        correlations_5p.append(corr_5p)
        pvalues_5p.append(pval_5p)
        correlations_3p.append(corr_3p)
        pvalues_3p.append(pval_3p)

    return correlations_5p, pvalues_5p, correlations_3p, pvalues_3p

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_correlations.py <input_matrix_file> <sample_info_file> <output_file>")
        sys.exit(1)

    input_matrix_file = sys.argv[1]
    sample_info_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load the matrix and sample info
    matrix_df = pd.read_csv(input_matrix_file, sep="\t")
    sample_info = load_sample_info(sample_info_file)

    # Calculate correlations
    correlations_5p, pvalues_5p, correlations_3p, pvalues_3p = calculate_correlations(matrix_df, sample_info)

    # Add new columns to the matrix
    matrix_df['Correlation_5p'] = correlations_5p
    matrix_df['Correlation_Pvalue_5p'] = pvalues_5p
    matrix_df['Correlation_3p'] = correlations_3p
    matrix_df['Correlation_Pvalue_3p'] = pvalues_3p

    # Save the final matrix with correlation results
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Matrix with correlation results has been created successfully in {output_file}")