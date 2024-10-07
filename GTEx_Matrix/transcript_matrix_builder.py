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
        pickle_path = os.path.join(info['file_path'], f"{srr}.ENST.result.pkl")
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as pickle_file:
                sample_dicts[srr] = pickle.load(pickle_file)
        else:
            print(f"Warning: Pickle file not found for {srr}: {pickle_path}")
    return sample_dicts

def read_refDB(refDB_file):
    chimeric_refs = []
    with open(refDB_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                ref = line.strip()[1:]  # Remove '>' at the beginning
                if 'x' in ref:
                    chimeric_refs.append(ref)
    return chimeric_refs

def process_chimeras(chimeric_refs, sample_dicts, sample_info):
    results = {}
    samples = list(sample_dicts.keys())
    
    for ref in chimeric_refs:
        row = [ref]
        for sample in samples:
            sample_dict = sample_dicts[sample]
            if ref in sample_dict:
                chi_count, p1_count, p2_count, ratio_5p, ratio_3p = map(float, sample_dict[ref].split(','))
                row.extend([chi_count, p1_count, p2_count, ratio_5p, ratio_3p])
            else:
                row.extend([0, 0, 0, 0, 0])
        results[ref] = row
    
    columns = ['Chimera']
    for sample in samples:
        columns.extend([f"{sample}_chi", f"{sample}_5p", f"{sample}_3p", 
                        f"{sample}_5_ratio", f"{sample}_3_ratio"])
    df = pd.DataFrame(list(results.values()), columns=columns)
    return df

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python build_matrix.py <sample_info_file> <refDB_file> <output_file>")
        sys.exit(1)
    
    sample_info_file = sys.argv[1]
    refDB_file = sys.argv[2]
    output_file = sys.argv[3]
    
    sample_info = load_sample_info(sample_info_file)
    sample_dicts = load_sample_dicts(sample_info)
    chimeric_refs = read_refDB(refDB_file)
    matrix_df = process_chimeras(chimeric_refs, sample_dicts, sample_info)
    
    matrix_df.to_csv(output_file, sep="\t", index=False, float_format='%.6f')
    print(f"Chimera matrix has been created successfully in {output_file}")