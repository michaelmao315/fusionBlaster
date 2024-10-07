import os
import sys
import glob

def combine_alignment_tsv(output_dir, fastq):
    # Get all readAlignment.tsv files for this fastq
    tsv_files = glob.glob(os.path.join(output_dir, f"{fastq}_chunk*.readAlignment.tsv"))
    
    # Output file
    output_file = os.path.join(output_dir, f"{fastq}.readAlignment.tsv")
    
    with open(output_file, 'w') as outfile:
        for tsv_file in tsv_files:
            with open(tsv_file, 'r') as infile:
                outfile.write(infile.read())
            
            # Remove the chunk file after processing
            os.remove(tsv_file)
            print(f"Processed and removed: {tsv_file}")
    
    print(f"Combined TSV saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python combine_alignment_tsv.py <output_dir> <fastq>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    fastq = sys.argv[2]
    combine_alignment_tsv(output_dir, fastq)