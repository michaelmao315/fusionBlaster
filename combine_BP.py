import os
import sys
import pickle
import glob

def combine_bp_results(output_dir, fastq):
    combined_results = {
        'chimera': {},
        'parental_5p': {},
        'parental_3p': {}
    }
    
    # Get all BP result files for this fastq
    bp_result_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*.genomic_bp_counts.pkl"))
    
    for file in bp_result_files:
        with open(file, 'rb') as f:
            chunk_results = pickle.load(f)
        
        # Combine the dictionaries for each category
        for category in ['chimera', 'parental_5p', 'parental_3p']:
            for reference, count in chunk_results[category].items():
                if reference in combined_results[category]:
                    combined_results[category][reference] += count
                else:
                    combined_results[category][reference] = count
    
      # Remove the chunk file after processing
        os.remove(file)

    # Save the combined BP results
    with open(os.path.join(output_dir, f"{fastq}.genomic_bp_counts.pkl"), 'wb') as f:
        pickle.dump(combined_results, f)
    
    print(f"Combined BP results saved as {fastq}.genomic_bp_counts.pkl")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_BP_results.py <output_dir> <fastq>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    fastq = sys.argv[2]
    combine_bp_results(output_dir, fastq)