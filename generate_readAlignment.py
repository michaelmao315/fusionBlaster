import pickle
import sys
import os
import pickle
import sys
import os

fastq = sys.argv[1]
output_dir = sys.argv[2]

# Load the readHashc dictionary
with open(os.path.join(output_dir, f"{fastq}.readHashc.pkl"), "rb") as f:
    readHashc = pickle.load(f)

# Create a new dictionary to store the alignment information
readAlignment = {}

with open(os.path.join(output_dir, f"{fastq}.readAlignment.tsv"), "w") as out2:
    # Iterate over each read in the readHashc dictionary
    for read in readHashc.keys():
        if readHashc[read]:
            # Find the longest read length by converting lengths to integers
            longest_length = max(int(length) for length in readHashc[read].keys())
            
            # Iterate over entries with the longest read length
            for ref in readHashc[read][longest_length]:  
                out2.write(f"{read}\t{ref}\n")
                
                # Add to readAlignment dictionary
                if read not in readAlignment:
                    readAlignment[read] = []
                readAlignment[read].append(ref)
                
# Convert lists to comma-separated strings
for read in readAlignment:
    readAlignment[read] = ",".join(readAlignment[read])

# Save the readAlignment dictionary as a pickle file
with open(os.path.join(output_dir, f"{fastq}.readAlignment.pkl"), "wb") as f:
    pickle.dump(readAlignment, f)

print(f"Processing complete. Results saved in {fastq}.readAlignment.tsv and {fastq}.readAlignment.pkl")