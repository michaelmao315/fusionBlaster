import pickle
import sys
import os

def generate_alignment_tsv(output_dir, readhashc_file, chunk_name):
    # Load the readHashc dictionary
    with open(readhashc_file, "rb") as f:
        readHashc = pickle.load(f)

    output_file = os.path.join(output_dir, f"{chunk_name}.readAlignment.tsv")
    
    with open(output_file, "w") as out:
        # Iterate over each read in the readHashc dictionary
        for read in readHashc.keys():
            if readHashc[read]:
                # Find the longest read length by converting lengths to integers
                longest_length = max(int(length) for length in readHashc[read].keys())
                
                # Iterate over entries with the longest read length
                for ref in readHashc[read][longest_length]:  
                    out.write(f"{read}\t{ref}\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_alignment_tsv.py <output_dir> <readhashc_file> <chunk_name>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    readhashc_file = sys.argv[2]
    chunk_name = sys.argv[3]
    generate_alignment_tsv(output_dir, readhashc_file, chunk_name)