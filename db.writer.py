#!/usr/bin/env python3
import sys
import subprocess
import os
import shutil

# Checking if enough arguments are passed
if len(sys.argv) < 3:
    print("Usage: script.py <input_file> <genome>")
    sys.exit(1)

# Reading command line arguments
input_file = sys.argv[1]
genome = sys.argv[2]

# Create a debug directory if it doesn't exist
debug_dir = "debug_files"
if not os.path.exists(debug_dir):
    os.makedirs(debug_dir)

ensts = set()
coords = []  # Store necessary data for main processing

with open(input_file, 'r') as file:
    for line in file:
        line = line.strip()
        coords.append(line.split('\t'))  # Store as parts, not raw lines
        for col in line.split('\t')[2:4]:  # Process only relevant columns
            ensts.update(part.split(',')[1] for part in col.split(':') if ',' in part)

ensts_str = '\n'.join(sorted(ensts))

def run_r_script_and_create_dict(script_path, genome, ensts_str):
    # Write ENSTs to a persistent file in the debug directory
    enst_file_path = os.path.join(debug_dir, "ensts_list.txt")
    with open(enst_file_path, 'w') as enst_file:
        enst_file.write(ensts_str)

    # Run the R script with the path to the ENST file
    result = subprocess.run(["Rscript", script_path, genome, enst_file_path], 
                          capture_output=True, text=True)

    # Process the output into a dictionary
    transcriptome_dict = {}
    output_lines = result.stdout.strip().split('\n')
    
    # Save R script output for debugging
    output_file_path = os.path.join(debug_dir, "r_script_output.txt")
    with open(output_file_path, 'w') as output_file:
        output_file.write(result.stdout)
        if result.stderr:
            output_file.write("\n=== STDERR ===\n")
            output_file.write(result.stderr)

    for line in output_lines:
        parts = line.split(maxsplit=1)  # Splitting the line into two parts: ID and sequence
        if len(parts) == 2:
            transcriptome_name_id, sequence = parts
            transcriptome_dict[transcriptome_name_id] = sequence
        else:
            print(f"Could not parse line: {line}")  # For lines that do not match the expected format

    return transcriptome_dict

transcriptome_script_path = os.path.expanduser("./fusionBlaster/transcriptome.R")
transcriptome = run_r_script_and_create_dict(transcriptome_script_path, genome, ensts_str)

unique_entries = {}  # To store unique parent and chimeric sequences

def update_or_add_entry(unique_entries, id, bp, seq):
    existing_key = next((key for key, value in unique_entries.items() if value == seq), None)
    if existing_key:
        existing_bps = existing_key.split("#")[1].split(",")
        new_bps = bp.split(",")
        all_bps = sorted(set(existing_bps + new_bps), key=int)
        new_name1 = f">{id}#{','.join(map(str, all_bps))}"
        if new_name1 != existing_key:
            unique_entries[new_name1] = seq
            del unique_entries[existing_key]
    else:
        unique_entries[f">{id}#{bp}"] = seq

for parts in coords:
    if len(parts) < 4:
        continue
    
    ENSG1, ENSG2 = parts[0], parts[1]
    
    for coord1 in parts[2].split(':'):
        if ',' not in coord1:
            continue
        bp1, enst1 = coord1.split(',')
        seq1 = transcriptome.get(enst1, "")
        if seq1:
            up_seq = seq1[:int(bp1)]
            update_or_add_entry(unique_entries, f"{ENSG1}:{enst1}", bp1, seq1)
            
            for coord2 in parts[3].split(':'):
                if ',' not in coord2:
                    continue
                bp2, enst2 = coord2.split(',')
                seq2 = transcriptome.get(enst2, "")
                if seq2:
                    update_or_add_entry(unique_entries, f"{ENSG2}:{enst2}", bp2, seq2)
                    chimeric_seq = up_seq + seq2[int(bp2) - 1:]
                    chimeric_name = f">{ENSG1}:{enst1}x{ENSG2}:{enst2},{bp2}#{bp1}"
                    unique_entries[chimeric_name] = chimeric_seq

# Save unique entries to a debug file
debug_output_path = os.path.join(debug_dir, "unique_entries.txt")
with open(debug_output_path, 'w') as debug_file:
    for entry_id, sequence in unique_entries.items():
        debug_file.write(f"{entry_id}\n{sequence}\n")

shutil.rmtree(debug_dir)

# Print results
for entry_id, sequence in unique_entries.items():
    print(entry_id)
    print(sequence)

print(">ENSG00000075624:ENST00000493945#1086")
print("CTCACCATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCCCCCGGGCCGTCTTCCCCTCCATCGTGGGGCGCCCCAGGCACCAGGGCGTGATGGTGGGCATGGGTCAGAAGGATTCCTATGTGGGCGACGAGGCCCAGAGCAAGAGAGGCATCCTCACCCTGAAGTACCCCATCGAGCACGGCATCGTCACCAACTGGGACGACATGGAGAAAATCTGGCACCACACCTTCTACAATGAGCTGCGTGTGGCTCCCGAGGAGCACCCCGTGCTGCTGACCGAGGCCCCCCTGAACCCCAAGGCCAACCGCGAGAAGATGACCCAGATCATGTTTGAGACCTTCAACACCCCAGCCATGTACGTTGCTATCCAGGCTGTGCTATCCCTGTACGCCTCTGGCCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACCCACACTGTGCCCATCTACGAGGGGTATGCCCTCCCCCATGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACGGCCGAGCGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAAGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGTGAGTGGAGACTGTCTCCCGGCTCTGCCTGACATGAGGGTTACCCCTCGGGGCTGTGCTGTGGAAGCTAAGTCCTGCCCTCATTTCCCTCTCAGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGCGGCACCACCATGTACCCTGGCATTGCCGACAGGATGCAGAAGGAGATCACTGCCCTGGCACCCAGCACAATGAAGATCAAGATCATTGCTCCTCCTGAGCGCAAGTACTCCGTGTGGATCGGCGGCTCCATCCTGGCCTCGCTGTCCACCTTCCAGCAGATGTGGATCAGCAAGCAGGAGTATGACGAGTCCGGCCCCTCCATCGTCCACCGCAAATGCTTCTAGGCGGACT")