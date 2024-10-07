import os
import sys
import pickle
import glob

def combine_geneHashcc(output_dir, fastq):
    combined_geneHashcc = {}
    
    # Get all geneHashcc files for this fastq
    geneHashcc_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk_*.geneHashcc.pkl"))
    
    for file in geneHashcc_files:
        with open(file, 'rb') as f:
            chunk_geneHashcc = pickle.load(f)
        
        # Combine the dictionaries
        for gene, count in chunk_geneHashcc.items():
            if gene in combined_geneHashcc:
                combined_geneHashcc[gene] += count
            else:
                combined_geneHashcc[gene] = count
        os.remove(file)
        
        
    # Save the combined geneHashcc
    with open(os.path.join(output_dir, f"{fastq}.geneHashcc.pkl"), 'wb') as f:
        pickle.dump(combined_geneHashcc, f)
    
    print(f"Combined geneHashcc saved as {fastq}.geneHashcc.pkl")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_geneHashcc_results.py <output_dir> <fastq>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    fastq = sys.argv[2]
    combine_geneHashcc(output_dir, fastq)