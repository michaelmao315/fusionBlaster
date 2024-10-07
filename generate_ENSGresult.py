import pickle
import sys
import os

fastq = sys.argv[1]
output_dir = sys.argv[2]

# Load the geneHashcc dictionary
with open(os.path.join(output_dir, f"{fastq}.geneHashcc.pkl"), "rb") as f:
    geneHashcc = pickle.load(f)

# Create the new dictionary for ENSG results
ENSG_result = {}

with open(os.path.join(output_dir, f"{fastq}.ENSG.result.tsv"), "w") as out_file:
    for key in sorted(geneHashcc.keys()):
        if 'x' in key:
            # This is a chimeric RNA
            g1, g2 = key.split('x')
            
            # Ensure g1 and g2 exist in geneHashcc
            if g1 not in geneHashcc:
                geneHashcc[g1] = 0
            if g2 not in geneHashcc:
                geneHashcc[g2] = 0
            
            # Calculate ratios with division by zero checks
            if geneHashcc[key] + geneHashcc[g1] > 0:
                r1 = geneHashcc[key] / (geneHashcc[key] + geneHashcc[g1])
            else:
                r1 = 0
            
            if geneHashcc[key] + geneHashcc[g2] > 0:
                r2 = geneHashcc[key] / (geneHashcc[key] + geneHashcc[g2])
            else:
                r2 = 0
            
            # Write to the output file
            out_line = f"{key}\t{geneHashcc[key]}\t{geneHashcc[g1]}\t{geneHashcc[g2]}\t{r1:.6f}\t{r2:.6f}\n"
            out_file.write(out_line)
            
            # Store in ENSG_result dictionary
            ENSG_result[key] = f"{geneHashcc[key]},{geneHashcc[g1]},{geneHashcc[g2]},{r1:.6f},{r2:.6f}"
            
# Save ENSG_result as a pickle file
with open(os.path.join(output_dir, f"{fastq}.ENSG.result.pkl"), "wb") as f:
    pickle.dump(ENSG_result, f)

print(f"Processing complete. Results saved in {fastq}.ENSG.result.tsv and {fastq}.ENSG.result.pkl")