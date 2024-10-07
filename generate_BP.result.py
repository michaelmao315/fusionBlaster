import sys
import os
import pickle

def load_results(result_file):
    with open(result_file, 'rb') as f:
        return pickle.load(f)

def process_input_file(input_file, results, output_file, pickle_output_file):
    chimera_counts = results["chimera"]
    parental_5p_counts = results['parental_5p']
    parental_3p_counts = results['parental_3p']
    
    output_dict = {}
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                # Combine the first 5 columns into one, separated by commas
                combined_key = ','.join(parts[:5])
                
                genomic_bp = f"{parts[3]}|{parts[4]}"
                chimera_count = chimera_counts.get(genomic_bp, 0)
                parental_5p_count = parental_5p_counts.get(genomic_bp, 0)
                parental_3p_count = parental_3p_counts.get(genomic_bp, 0)
                
                # Calculate ratios
                five_p_ratio = chimera_count / (chimera_count + parental_5p_count) if chimera_count + parental_5p_count > 0 else 0
                three_p_ratio = chimera_count / (chimera_count + parental_3p_count) if chimera_count + parental_3p_count > 0 else 0
                
                # Construct new line with combined key and additional columns
                new_line = f"{combined_key}\t{chimera_count}\t{parental_5p_count}\t{parental_3p_count}\t{five_p_ratio:.6f}\t{three_p_ratio:.6f}\n"
                outfile.write(new_line)
                
                # Add to dictionary
                output_dict[combined_key] = f"{chimera_count},{parental_5p_count},{parental_3p_count},{five_p_ratio:.6f},{three_p_ratio:.6f}"
            else:
                outfile.write(line)  # Write unchanged if not enough columns
    
    # Save the dictionary as a pickle file
    with open(pickle_output_file, 'wb') as pickle_file:
        pickle.dump(output_dict, pickle_file)

def main(fastq, output_dir, input_file):
    result_file = os.path.join(output_dir, f"{fastq}.genomic_bp_counts.pkl")
    output_file = os.path.join(output_dir, f"{fastq}.GenoBP.result.tsv")
    pickle_output_file = os.path.join(output_dir, f"{fastq}.GenoBP.result.pkl")
    
    results = load_results(result_file)
    process_input_file(input_file, results, output_file, pickle_output_file)
    print(f"Processing complete. TSV results saved in {output_file}")
    print(f"Dictionary results saved in {pickle_output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python /project/hlilab/software/fusionBlaster/bp_merge.py <fastq> <output_dir> <input_file>")
        sys.exit(1)
    
    fastq = sys.argv[1]
    output_dir = sys.argv[2]
    input_file = sys.argv[3]
    
    main(fastq, output_dir, input_file)