import pickle
import sys
import os

fastq = sys.argv[1]
output_dir = sys.argv[2]

with open(os.path.join(output_dir, f"{fastq}.readHashcc.pkl"), "rb") as f:
    readHashcc = pickle.load(f)

# Create the new dictionary
ENST_result = {}

with open(os.path.join(output_dir, f"{fastq}.ENST.result.tsv"), "w") as out1:
    for col23 in sorted(readHashcc.keys()):
        if 'x' in col23:
            g1, g2 = col23.split('x')
            cpg2 = g2.split(',')[-1]
            three, five = cpg2.split('#')
            g2 = g2.rsplit(',', 1)[0] + '#' + three
            g1 = g1 + '#' + five
            
            if g1 not in readHashcc:
                readHashcc[g1] = 0
            if g2 not in readHashcc:
                readHashcc[g2] = 0
            
            # Avoid division by zero
            if readHashcc[col23] + readHashcc[g1] > 0:
                r5 = readHashcc[col23] / (readHashcc[col23] + readHashcc[g1])
            else:
                r5 = 0
            
            if readHashcc[col23] + readHashcc[g2] > 0:
                r3 = readHashcc[col23] / (readHashcc[col23] + readHashcc[g2])
            else:
                r3 = 0
            
            # Write to result1.tsv
            out1.write(f"{col23}\t{readHashcc[col23]}\t{readHashcc[g1]}\t{readHashcc[g2]}\t{r5:.6f}\t{r3:.6f}\n")
            
            # Store in ENST_result dictionary
            ENST_result[col23] = f"{readHashcc[col23]},{readHashcc[g1]},{readHashcc[g2]},{r5:.6f},{r3:.6f}"
            
# Save ENST_result as a pickle file
with open(os.path.join(output_dir, f"{fastq}.ENST.result.pkl"), "wb") as f:
    pickle.dump(ENST_result, f)