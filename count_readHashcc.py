import pickle
import sys
import os

def main(output_dir, readhashc_file, chunk_name):
    # Load the readHashc dictionary
    with open(readhashc_file, "rb") as f:
        readHashc = pickle.load(f)

    # Dictionary to store read counts
    readHashcc = {}

    # Process each key in readHashc
    for read in readHashc.keys():
        # Ensure there are read lengths to process
        if readHashc[read]:
            # Find the longest read length
            longest_length = max(int(length) for length in readHashc[read].keys())
            
            # Process entries for the longest read length
            for entry in readHashc[read][longest_length]:
                readHashcc.setdefault(entry, 0)
                readHashcc[entry] += 1

    # Save the updated readHashcc dictionary
    with open(os.path.join(output_dir, f"{chunk_name}.readHashcc.pkl"), "wb") as f:
        pickle.dump(readHashcc, f)

    print(f"Processing complete for {chunk_name}. Results saved in {chunk_name}.readHashcc.pkl")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <output_dir> <readhashc_file> <chunk_name>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])