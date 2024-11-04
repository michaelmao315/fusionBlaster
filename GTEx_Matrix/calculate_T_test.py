import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def load_sample_info(sample_info_file):
    sample_info = {}
    with open(sample_info_file, 'r') as f:
        for line in f:
            srr, group, _ = line.strip().split('\t')
            sample_info[srr] = group
    return sample_info

def add_T_pvalues(matrix_df, group1, group2, ratio_type):
    pvalues = []
    mean_differences = []
    fold_changes = []
    group1_means = []
    group2_means = []

    for _, row in matrix_df.iterrows():
        group1_values = [float(row[f"{sample}_{ratio_type}_ratio"]) for sample in group1]
        group2_values = [float(row[f"{sample}_{ratio_type}_ratio"]) for sample in group2]
        
        _, pvalue = ttest_ind(group1_values, group2_values)
        group1_mean = np.mean(group1_values)
        group2_mean = np.mean(group2_values)
        mean_difference = group1_mean - group2_mean

        if group1_mean == 0 and group2_mean == 0:
            fold_change = 1
        else:
            fold_change = group1_mean / group2_mean if group2_mean != 0 else np.inf


        pvalues.append(pvalue if not np.isnan(pvalue) else 1)
        mean_differences.append(mean_difference)
        fold_changes.append(fold_change)
        group1_means.append(group1_mean)
        group2_means.append(group2_mean)

    corrected_pvalues = multipletests(pvalues, method='fdr_bh')[1]
    
    log2_fold_changes = [np.log2(fc) if fc > 0 else 0 for fc in fold_changes]
    neg_log10_pvalues = [-np.log10(p) if p > 0 else 10 for p in pvalues]
    neg_log10_corrected_pvalues = [-np.log10(p) if p > 0 else 10 for p in corrected_pvalues]
    
    return pvalues, corrected_pvalues, mean_differences, fold_changes, log2_fold_changes, neg_log10_pvalues, neg_log10_corrected_pvalues, group1_means, group2_means

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python calculate_ttest.py <matrix_file> <sample_info_file> <output_file> <group1_name> <group2_name>")
        sys.exit(1)

    matrix_file = sys.argv[1]
    sample_info_file = sys.argv[2]
    output_file = sys.argv[3]
    group1_name = sys.argv[4]
    group2_name = sys.argv[5]

    # Load the matrix and sample info
    matrix_df = pd.read_csv(matrix_file, sep="\t")
    sample_info = load_sample_info(sample_info_file)
    chimera_column = matrix_df.columns[0]
    matrix_df = matrix_df.drop_duplicates(subset=[chimera_column], keep='first')

    # Determine the groups based on sample info and user input
    group1 = [srr for srr, group in sample_info.items() if group == group1_name]
    group2 = [srr for srr, group in sample_info.items() if group == group2_name]

    if not group1 or not group2:
        print(f"Error: One or both groups ({group1_name}, {group2_name}) not found in sample info.")
        sys.exit(1)

    # Calculate T-test for 5' ratio
    pvalues_5p, corrected_pvalues_5p, mean_differences_5p, fold_changes_5p, log2_fold_changes_5p, neg_log10_pvalues_5p, neg_log10_corrected_pvalues_5p, group1_means_5p, group2_means_5p = add_T_pvalues(matrix_df, group1, group2, '5')

    # Calculate T-test for 3' ratio
    pvalues_3p, corrected_pvalues_3p, mean_differences_3p, fold_changes_3p, log2_fold_changes_3p, neg_log10_pvalues_3p, neg_log10_corrected_pvalues_3p, group1_means_3p, group2_means_3p = add_T_pvalues(matrix_df, group1, group2, '3')

    # Prepare the new columns
    new_columns = {
        f'T_Pvalue_5p': pvalues_5p,
        f'T_Corrected_Pvalue_5p': corrected_pvalues_5p,
        f'MeanDifference_{group1_name}_vs_{group2_name}_5p': mean_differences_5p,
        f'FoldChange_{group1_name}_vs_{group2_name}_5p': fold_changes_5p,
        f'Log2FoldChange_{group1_name}_vs_{group2_name}_5p': log2_fold_changes_5p,
        f'NegLog10Pvalue_5p': neg_log10_pvalues_5p,
        f'NegLog10CorrectedPvalue_5p': neg_log10_corrected_pvalues_5p,
        f'Mean_{group1_name}_5p': group1_means_5p,
        f'Mean_{group2_name}_5p': group2_means_5p,
        f'T_Pvalue_3p': pvalues_3p,
        f'T_Corrected_Pvalue_3p': corrected_pvalues_3p,
        f'MeanDifference_{group1_name}_vs_{group2_name}_3p': mean_differences_3p,
        f'FoldChange_{group1_name}_vs_{group2_name}_3p': fold_changes_3p,
        f'Log2FoldChange_{group1_name}_vs_{group2_name}_3p': log2_fold_changes_3p,
        f'NegLog10Pvalue_3p': neg_log10_pvalues_3p,
        f'NegLog10CorrectedPvalue_3p': neg_log10_corrected_pvalues_3p,
        f'Mean_{group1_name}_3p': group1_means_3p,
        f'Mean_{group2_name}_3p': group2_means_3p
    }

    # Add new columns to the existing DataFrame
    insert_position = 1
    for col_name, col_data in new_columns.items():
        matrix_df.insert(insert_position, col_name, col_data)
        insert_position += 1

    # Save the final matrix with T-test results
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Matrix with T-test results has been created successfully in {output_file}")
