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
        pickle_path = os.path.join(info['file_path'], f"{srr}.GenoBP.result.pkl")
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as pickle_file:
                sample_dicts[srr] = pickle.load(pickle_file)
        else:
            print(f"Warning: Pickle file not found for {srr}: {pickle_path}")
    return sample_dicts

def process_chimeras(chimera_file, sample_dicts, sample_info):
    results = {}
    samples = list(sample_dicts.keys())
    
    with open(chimera_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 5:
                chimera_ref = ','.join(parts[:5])
                if chimera_ref not in results:
                    row = [chimera_ref]
                    for sample in samples:
                        sample_dict = sample_dicts[sample]
                        if chimera_ref in sample_dict:
                            chi_count, p5_count, p3_count, ratio_5p, ratio_3p = map(float, sample_dict[chimera_ref].split(','))
                            row.extend([chi_count, p5_count, p3_count, ratio_5p, ratio_3p])
                        else:
                            row.extend([0, 0, 0, 0, 0])
                    results[chimera_ref] = row
    
    columns = ['Chimera']
    for sample in samples:
        columns.extend([f"{sample}_chi", f"{sample}_5p", f"{sample}_3p", 
                        f"{sample}_5_ratio", f"{sample}_3_ratio"])
    df = pd.DataFrame(list(results.values()), columns=columns)
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