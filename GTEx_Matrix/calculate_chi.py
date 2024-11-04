import sys
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

def load_sample_info(sample_info_file):
    sample_info = {}
    with open(sample_info_file, 'r') as f:
        for line in f:
            srr, group, _ = line.strip().split('\t')
            sample_info[srr] = group
    return sample_info

def add_ChiSquare_pvalues(matrix_df, group_yes, group_no):
    pvalues = []

    for _, row in matrix_df.iterrows():
        groupYes = [1 if float(row[f"{sample}_chi"]) > 0 else 0 for sample in group_yes]
        groupNo = [1 if float(row[f"{sample}_chi"]) > 0 else 0 for sample in group_no]

        contingency_table = [
            [sum(groupYes), len(groupYes) - sum(groupYes)],
            [sum(groupNo), len(groupNo) - sum(groupNo)]
        ]

        pvalue = calculate_pvalue(contingency_table)
        pvalues.append(pvalue)

    corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1]
    neg_log10_pvalues = [-np.log10(p) if p > 0 else 10 for p in pvalues]
    neg_log10_corrected_pvalues = [-np.log10(p) if p > 0 else 10 for p in corrected_pvalues]

    return (pvalues, corrected_pvalues, neg_log10_pvalues, neg_log10_corrected_pvalues)

def calculate_pvalue(contingency_table):
    contingency_table = np.array(contingency_table)
    if np.any(np.all(contingency_table == 0, axis=0)):
        return 1.0  # Assign a p-value of 1 if any column has all zeros
    elif np.all(contingency_table == contingency_table[0]):
        return 1.0  # Assign a p-value of 1 if all rows are the same
    else:
        _, pvalue, _, _ = chi2_contingency(contingency_table, correction=False)
        return pvalue

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python calculate_chi_square.py <input_matrix_file> <sample_info_file> <output_file> <group1_name> <group2_name>")
        sys.exit(1)

    input_matrix_file = sys.argv[1]
    sample_info_file = sys.argv[2]
    output_file = sys.argv[3]
    group1_name = sys.argv[4]
    group2_name = sys.argv[5]

    # Load the matrix and sample info
    matrix_df = pd.read_csv(input_matrix_file, sep="\t")
    sample_info = load_sample_info(sample_info_file)

    # Determine the groups based on sample info
    group_yes = [srr for srr, group in sample_info.items() if group == group1_name]
    group_no = [srr for srr, group in sample_info.items() if group == group2_name]

    # Calculate Chi-square test results
    pvalues, corrected_pvalues, neg_log10_pvalues, neg_log10_corrected_pvalues = add_ChiSquare_pvalues(matrix_df, group_yes, group_no)

    # Add the new columns to the DataFrame
    new_columns = {
        'ChiSquare_Pvalue': pvalues,
        'ChiSquare_Corrected_Pvalue': corrected_pvalues,
        'ChiSquare_NegLog10Pvalue': neg_log10_pvalues,
        'ChiSquare_NegLog10CorrectedPvalue': neg_log10_corrected_pvalues
    }

    # Insert the new columns after the first column
    insert_position = 1
    for col_name, col_data in new_columns.items():
        matrix_df.insert(insert_position, col_name, col_data)
        insert_position += 1  # Increment the position for the next column

    # Save the final matrix with Chi-square test results
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Matrix with Chi-square test results has been created successfully in {output_file}")