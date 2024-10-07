import sys
import os
import pickle
import pandas as pd

def load_sample_info(sample_info_file):
    sample_info = {}
    with open(sample_info_file, 'r') as f:
        for line in f:
            srr, group, file_path = line.strip().split('\t')
            sample_info[srr] = {'group': group, 'file_path': file_path}
    return sample_info

def load_sample_dicts(sample_info):
    sample_dicts = {}
    for srr, info in sample_info.items():
        pickle_path = os.path.join(info['file_path'], f"{srr}.geneHashcc.pkl")
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as pickle_file:
                sample_dicts[srr] = pickle.load(pickle_file)
        else:
            print(f"Warning: Pickle file not found for {srr}: {pickle_path}")
    return sample_dicts

def process_chimeras(chimera_file, sample_dicts, sample_info):
    results = []
    samples = list(sample_dicts.keys())
    
    with open(chimera_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                ensg1, ensg2 = parts[1], parts[2]
                chimera = f"{ensg1}x{ensg2}"
                row = [chimera]
                for sample in samples:
                    sample_dict = sample_dicts[sample]
                    chi_count = sample_dict.get(chimera, 0)
                    p1_count = sample_dict.get(ensg1, 0)
                    p2_count = sample_dict.get(ensg2, 0)
                    
                    ratio_5p = chi_count / (p1_count+chi_count) if chi_count != 0 else 0
                    ratio_3p = chi_count / (p2_count+chi_count) if chi_count != 0 else 0
                    row.extend([chi_count, p1_count, p2_count, ratio_5p, ratio_3p])
                results.append(row)
    
    columns = ['Chimera']
    for sample in samples:
        columns.extend([f"{sample}_chi", f"{sample}_p1", f"{sample}_p2", 
                        f"{sample}_5_ratio", f"{sample}_3_ratio"])
    df = pd.DataFrame(results, columns=columns)
    return df

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python build_matrix.py <sample_info_file> <chimera_file> <output_file>")
        sys.exit(1)
    
    sample_info_file = sys.argv[1]
    chimera_file = sys.argv[2]
    output_file = sys.argv[3]
    
    sample_info = load_sample_info(sample_info_file)
    sample_dicts = load_sample_dicts(sample_info)
    matrix_df = process_chimeras(chimera_file, sample_dicts, sample_info)
    
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Chimera matrix has been created successfully in {output_file}")