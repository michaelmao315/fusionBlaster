import pickle
import sys
import os
import re

def preprocess_breakpoint_info(breakpoint_file):
    breakpoint_info = {}
    gene_to_chimera = {'5p': {}, '3p': {}}
    chimera_lookup = {}

    with open(breakpoint_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            gene_pair = f"{parts[1]}x{parts[2]}"
            genomic_bp = f"{parts[3]}|{parts[4]}"

            breakpoint_info[genomic_bp] = {
                '5p': {trx.split(',')[1]: trx.split(',')[0] for trx in parts[7].split(':') if ',' in trx},
                '3p': {trx.split(',')[1]: trx.split(',')[0] for trx in parts[8].split(':') if ',' in trx}
            }

            gene_to_chimera['5p'].setdefault(parts[1], set()).add(genomic_bp)
            gene_to_chimera['3p'].setdefault(parts[2], set()).add(genomic_bp)

            for trx_5p, bp_5p in breakpoint_info[genomic_bp]['5p'].items():
                for trx_3p, bp_3p in breakpoint_info[genomic_bp]['3p'].items():
                    chimera_lookup[f"{gene_pair}_{trx_5p}_{trx_3p}_{bp_5p}_{bp_3p}"] = genomic_bp

    return breakpoint_info, gene_to_chimera, chimera_lookup

def get_gene_and_bp(ref):
    match = re.match(r'(ENSG\d+):(ENST\d+)(?:x(ENSG\d+):(ENST\d+),(\d+))?#(\d+)', ref)
    if match:
        if match.group(3):  # Chimeric case
            return f"{match.group(1)}x{match.group(3)}", match.group(2), match.group(4), match.group(6), match.group(5)
        else:  # Parental case
            return match.group(1), match.group(2), None, match.group(6), None
    return None, None, None, None, None


def process_read(read_data, breakpoint_info, gene_to_chimera, chimera_lookup):
    chimera_bps = set()
    parental_5p_bps = set()
    parental_3p_bps = set()

    if read_data:
        # Find the longest read length by using a generator expression for converting to integers
        longest_length = max(int(length) for length in read_data.keys())
        
        # Get the references for the longest read length
        references = read_data[longest_length]

        for ref in references:
            gene_key, trx_5p, trx_3p, bp_5p, bp_3p = get_gene_and_bp(ref)

            if trx_3p:  # Chimeric case
                lookup_key = f"{gene_key}_{trx_5p}_{trx_3p}_{bp_5p}_{bp_3p}"
                if lookup_key in chimera_lookup:
                    chimera_bps.add(chimera_lookup[lookup_key])
            else:  # Parental case
                if gene_key in gene_to_chimera['5p']:
                    for genomic_bp in gene_to_chimera['5p'][gene_key]:
                        if trx_5p in breakpoint_info[genomic_bp]['5p'] and bp_5p == breakpoint_info[genomic_bp]['5p'][trx_5p]:
                            parental_5p_bps.add(genomic_bp)
                if gene_key in gene_to_chimera['3p']:
                    for genomic_bp in gene_to_chimera['3p'][gene_key]:
                        if trx_5p in breakpoint_info[genomic_bp]['3p'] and bp_5p == breakpoint_info[genomic_bp]['3p'][trx_5p]:
                            parental_3p_bps.add(genomic_bp)

    return chimera_bps, parental_5p_bps, parental_3p_bps

def main(output_dir, breakpoint_file, readhashc_file, chunk_name):
    breakpoint_info, gene_to_chimera, chimera_lookup = preprocess_breakpoint_info(breakpoint_file)

    with open(readhashc_file, "rb") as f:
        readHashc = pickle.load(f)

    chimera_counts = {}
    parental_5p_counts = {}
    parental_3p_counts = {}

    for read in readHashc:
        c, p5, p3 = process_read(readHashc[read], breakpoint_info, gene_to_chimera, chimera_lookup)
        for genomic_bp in c:
            chimera_counts[genomic_bp] = chimera_counts.get(genomic_bp, 0) + 1
        for genomic_bp in p5:
            parental_5p_counts[genomic_bp] = parental_5p_counts.get(genomic_bp, 0) + 1
        for genomic_bp in p3:
            parental_3p_counts[genomic_bp] = parental_3p_counts.get(genomic_bp, 0) + 1

    results = {
        'chimera': chimera_counts,
        'parental_5p': parental_5p_counts,
        'parental_3p': parental_3p_counts
    }

    with open(os.path.join(output_dir, f"{chunk_name}.genomic_bp_counts.pkl"), "wb") as f:
        pickle.dump(results, f)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <output_dir> <breakpoint_file> <readhashc_file> <chunk_name>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
