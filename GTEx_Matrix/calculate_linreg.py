import sys
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def load_sample_info(sample_info_file):
    sample_info = {}
    with open(sample_info_file, 'r') as f:
        for line in f:
            srr, value, _ = line.strip().split('\t')
            sample_info[srr] = float(value)  # Convert to float for numerical analysis
    return sample_info

def calculate_linear_regression(matrix_df, sample_info):
    r_squared_5p = []
    pvalues_5p = []
    slopes_5p = []
    r_squared_3p = []
    pvalues_3p = []
    slopes_3p = []

    for _, row in matrix_df.iterrows():
        x = []  # age or BMI
        y_5p = []  # 5' ratios
        y_3p = []  # 3' ratios

        for sample in sample_info.keys():
            if f"{sample}_5_ratio" in row and f"{sample}_3_ratio" in row:
                x.append(sample_info[sample])
                y_5p.append(float(row[f"{sample}_5_ratio"]))
                y_3p.append(float(row[f"{sample}_3_ratio"]))

        # Linear regression for 5' ratio
        slope_5p, intercept_5p, r_value_5p, p_value_5p, std_err_5p = stats.linregress(x, y_5p)
        r_squared_5p.append(r_value_5p**2)
        pvalues_5p.append(p_value_5p)
        slopes_5p.append(slope_5p)

        # Linear regression for 3' ratio
        slope_3p, intercept_3p, r_value_3p, p_value_3p, std_err_3p = stats.linregress(x, y_3p)
        r_squared_3p.append(r_value_3p**2)
        pvalues_3p.append(p_value_3p)
        slopes_3p.append(slope_3p)

    # Calculate adjusted p-values
    _, adjusted_pvalues_5p, _, _ = multipletests(pvalues_5p, method='fdr_bh')
    _, adjusted_pvalues_3p, _, _ = multipletests(pvalues_3p, method='fdr_bh')

    return r_squared_5p, pvalues_5p, adjusted_pvalues_5p, slopes_5p, r_squared_3p, pvalues_3p, adjusted_pvalues_3p, slopes_3p

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_linear_regression.py <input_matrix_file> <sample_info_file> <output_file>")
        sys.exit(1)

    input_matrix_file = sys.argv[1]
    sample_info_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load the matrix and sample info
    matrix_df = pd.read_csv(input_matrix_file, sep="\t")
    sample_info = load_sample_info(sample_info_file)

    # Calculate linear regression
    r_squared_5p, pvalues_5p, adjusted_pvalues_5p, slopes_5p, r_squared_3p, pvalues_3p, adjusted_pvalues_3p, slopes_3p = calculate_linear_regression(matrix_df, sample_info)

    # Prepare the new columns
    new_columns = [
        ('Slope_3p', slopes_3p),
        ('P_value_3p', pvalues_3p),
        ('Adjusted_P_value_3p', adjusted_pvalues_3p),
        ('R_squared_3p', r_squared_3p),
        ('Slope_5p', slopes_5p),
        ('P_value_5p', pvalues_5p),
        ('Adjusted_P_value_5p', adjusted_pvalues_5p),
        ('R_squared_5p', r_squared_5p)
    ]

    # Insert the new columns after the first column
    for i, (col_name, col_data) in enumerate(reversed(new_columns), 1):
        matrix_df.insert(1, col_name, col_data)

    # Save the final matrix with linear regression results
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Matrix with linear regression results has been created successfully in {output_file}")