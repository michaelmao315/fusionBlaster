import pickle
import sys
import os
import re

def get_gene_key(key):
    ensgs = re.findall(r'(ENSG\d+)', key)
    if len(ensgs) == 1:
        return ensgs[0]  # Parental gene case
    elif len(ensgs) == 2:
        return f"{ensgs[0]}x{ensgs[1]}"  # Chimeric RNA case
    else:
        return None  # Invalid case

def main(output_dir, readhashc_file, chunk_name):
    # Load the existing readHashc dictionary
    with open(readhashc_file, "rb") as f:
        readHashc = pickle.load(f)

    # Create a new dictionary to store gene and chimera counts
    geneHashcc = {}

    # Iterate through the reads
    for read in readHashc.keys():
        # Set to keep track of unique gene keys for this read
        read_gene_set = set()
        
        # Ensure read lengths are treated as integers
        if readHashc[read]:
            try:
                # Convert read lengths to integers
                longest_length = max(int(length) for length in readHashc[read].keys())
            except ValueError:
                # Handle case where conversion fails (e.g., non-numeric keys)
                print(f"Warning: Non-numeric read length encountered for read '{read}'")
                continue
            
            # Process only the longest read length
            for key in readHashc[read][longest_length]:
                # Extract gene key from the key
                gene_key = get_gene_key(key)
                
                if gene_key:
                    # Add gene key to the set for this read
                    read_gene_set.add(gene_key)

        # Count each unique gene or chimera for this read
        for gene_key in read_gene_set:
            geneHashcc.setdefault(gene_key, 0)
            geneHashcc[gene_key] += 1

    # Save the new geneHashcc dictionary
    with open(os.path.join(output_dir, f"{chunk_name}.geneHashcc.pkl"), "wb") as f:
        pickle.dump(geneHashcc, f)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <output_dir> <readhashc_file> <chunk_name>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])