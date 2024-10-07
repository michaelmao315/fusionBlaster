import os
import sys
import pickle
import glob

def combine_readHashcc(output_dir, fastq):
    combined_readHashcc = {}
    
    # Get all readHashcc files for this fastq
    readHashcc_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*.readHashcc.pkl"))
    
    for file in readHashcc_files:
        with open(file, 'rb') as f:
            chunk_readHashcc = pickle.load(f)
        
        # Combine the dictionaries
        for read, count in chunk_readHashcc.items():
            if read in combined_readHashcc:
                combined_readHashcc[read] += count
            else:
                combined_readHashcc[read] = count
                
        os.remove(file)
    
    # Save the combined readHashcc
    with open(os.path.join(output_dir, f"{fastq}.readHashcc.pkl"), 'wb') as f:
        pickle.dump(combined_readHashcc, f)
    
    print(f"Combined readHashcc saved as {fastq}.readHashcc.pkl")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_readHashcc_results.py <output_dir> <fastq>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    fastq = sys.argv[2]
    combine_readHashcc(output_dir, fastq)