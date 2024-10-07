import pickle
import sys
import os

def generate_alignment_dict(output_dir, readhashc_file, chunk_name):
    print("Warning: This script may consume a large amount of memory. Consider using the TSV output instead if memory issues occur.")
    
    # Load the readHashc dictionary
    with open(readhashc_file, "rb") as f:
        readHashc = pickle.load(f)

    # Create a new dictionary to store the alignment information
    readAlignment = {}

    # Iterate over each read in the readHashc dictionary
    for read in readHashc.keys():
        if readHashc[read]:
            # Find the longest read length by converting lengths to integers
            longest_length = max(int(length) for length in readHashc[read].keys())
            
            # Iterate over entries with the longest read length
            for ref in readHashc[read][longest_length]:  
                if read not in readAlignment:
                    readAlignment[read] = []
                readAlignment[read].append(ref)

    # Convert lists to comma-separated strings
    for read in readAlignment:
        readAlignment[read] = ",".join(readAlignment[read])

    # Save the readAlignment dictionary as a pickle file
    output_file = os.path.join(output_dir, f"{chunk_name}.readAlignment.pkl")
    with open(output_file, "wb") as f:
        pickle.dump(readAlignment, f)

    print(f"Dictionary output generated: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_alignment_dict.py <output_dir> <readhashc_file> <chunk_name>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    readhashc_file = sys.argv[2]
    chunk_name = sys.argv[3]
    generate_alignment_dict(output_dir, readhashc_file, chunk_name)