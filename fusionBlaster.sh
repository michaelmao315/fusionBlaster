#!/bin/bash

# Get the directory where the script is located, resolving symlinks
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
FUSION_BLASTER_DIR="$(dirname "$SCRIPT_DIR")"

in=$1
out=$2
fastq=$3
refDB=$4
chimerabp=$5
cores=$6


if [[ $in == "" || $out == "" || $fastq == "" || $refDB == "" ]]; then 
echo "Proper Usage:"
echo "$0 /full/path/to/directory/of/fastqs /full/path/to/output/directory uniqFASTQid full/path/to/refDB.fa /full/path/to/chimeras.input.tsv cores"
echo '=================================================================='
echo '/full/path/to/directory/of/fastqs should be full path to a central directory where all fastq files to be analyzed are located'
echo '/full/path/to/output/directory should be single output directory where results from all fastqs anlayzed by same chimeras.input.tsv are located'
echo 'uniqFASTQid should be used from fastq files as [uniqFASTQid]_R1.fq.gz or [uniqFASTQid]_2.fastq.gz'
echo '/full/path/to/refDB.fa should be full path to the refDB.fa file created as STDOUT by the refDBwriter.sh'
echo '/full/path/to/chimeras.input.tsv is the full path to the original input file with chimera names, ENSGs, and breakpoint info'
echo 'cores is the number of CPU cores to use'
exit 1
fi

fq1=$(ls $in | grep "^$fastq"_ | grep -f <( echo fq ; echo fastq) |  sort -u | head -1)
fq2=$(ls $in | grep "^$fastq"_ | grep -f <( echo fq ; echo fastq) |  sort -u | tail -1)
if [[ $fq1 == "" || $fq2 == "" || $fq1 != *".gz" || $fq2 != *".gz" ]]; then 
echo 'fastq files were either not found in the specified input directory or are not named correctly'
echo 'fastq files must be paired end fastq files name like:'
echo 'uniq_1.fastq.gz uniq_2.fastq.gz or uniq_1.fq.gz uniq_2.fq.gz'
echo 'uniq_R1.fastq.gz uniq_R2.fastq.gz or uniq_R1.fq.gz uniq_R2.fq.gz'
exit 1
fi

echo ${in}/$fq1 ${in}/$fq2

readLength=$(for read in $(gzcat ${in}/$fq1 | tail -n +2 | awk 'NR%4==1' | head -2000); do echo $read | wc -m ; done | awk '{ sum += $1-1 } END { if (NR > 0) print sum / NR }' | cut -d. -f1)

echo $readLength
echo $SCRIPT_DIR
echo $out

if [ ! -d "$out" ]; then
  mkdir -p "$out"
fi

#convert fastqs to a fasta with read_1 and read_2 pasted together
gzcat ${in}/$fq1 | paste - - - - | cut -f1,2 | tr "@""\t" ">""\n" > ${out}/${fastq}_tmp1
paste <(gzcat ${in}/$fq2 | paste - - - - | cut -f1 | tr "@" ">") <(gzcat ${in}/$fq2 | paste - - - - | cut -f2 | tr ACTG TGAC | rev) | tr "\t" "\n" > ${out}/${fastq}_tmp2
paste ${out}/${fastq}_tmp1 ${out}/${fastq}_tmp2 | awk -v IFS="\t" -v OFS="" '{if ($1~">") print $1; else print $1,$2 }' > ${out}/${fastq}_query.fa
rm -rf ${out}/${fastq}_tmp1 ${out}/${fastq}_tmp2

echo "Finished everything before pBlat"


#run pBLAT and filter alignments

split -l 100000    ${out}/${fastq}_query.fa ${out}/${fastq}_chunk_
rm ${out}/${fastq}_query.fa

eval "$(conda shell.bash hook)"

# Main processing loop
for chunk in ${out}/${fastq}_chunk_*; do
    chunk_prefix=$(basename "$chunk" .fa)
    
    echo "Processing chunk: $chunk_prefix"

    conda activate pblat_env

    # Run pBLAT on the chunk, directing output to ${out} directory
    pblat -fastMap -threads=$cores -out=blast8 $refDB $chunk ${out}/${chunk_prefix}_pblat.out

    conda deactivate
    # Filter, directing blast output to ${out} directory
    awk -v OFS='\t' -v L=$readLength '{ if ($3>98 && $4>L-3) print $1,$2,$7,$8,$9,$10 }' ${out}/${chunk_prefix}_pblat.out > ${out}/${chunk_prefix}_blast.tsv
    
    # Remove the chunk and temporary pblat.out
    rm $chunk ${out}/${chunk_prefix}_pblat.out
    
    # Ensure every chunk of blast.tsv is below 10000000*cores lines
    lines=$(wc -l < ${out}/${chunk_prefix}_blast.tsv)
    if [ $lines -gt $((10000000 * cores)) ]; then
        python "${SCRIPT_DIR}/split_file.py" ${out}/${chunk_prefix}_blast.tsv ${out}/${chunk_prefix}_blast_division_ $((10000000 * cores))
        rm ${out}/${chunk_prefix}_blast.tsv
    else
        mv ${out}/${chunk_prefix}_blast.tsv ${out}/${chunk_prefix}_blast_division_0001.tsv
    fi

    # Run fusionBlaster_1.0.py on all blast TSV divisions
    for file in ${out}/${chunk_prefix}_blast_division_*.tsv; do
        division=$(basename "$file" .tsv | rev | cut -d '_' -f 1 | rev)
        python "${SCRIPT_DIR}/fusionBlaster.2.0.py" $file "${chunk_prefix}_${division}" "$out"
    done
    
    # Remove temporary files
    rm ${out}/${chunk_prefix}_blast_division_*.tsv
    
    echo "Finished processing chunk: $chunk_prefix"
done

echo "All chunks processed"

echo "Finished pBlat Loop"

for readhashc_file in ${out}/${fastq}_chunk*readHashc.pkl; do
    if [ -f "$readhashc_file" ]; then
        chunk_name=$(basename "$readhashc_file" .readHashc.pkl)
        
        # Run count_BP.py for each chunk
        python "${SCRIPT_DIR}/count_BP.py" "$out" "$chimerabp" "$readhashc_file" "$chunk_name"

        # Run count_readHashcc.py for each chunk
        python "${SCRIPT_DIR}/count_readHashcc.py" "$out" "$readhashc_file" "$chunk_name"

        # Run count_geneHashcc.py for each chunk
        python "${SCRIPT_DIR}/count_geneHashcc.py" "$out" "$readhashc_file" "$chunk_name"

        # Generate alignment TSV for each chunk
        python "${SCRIPT_DIR}/generate_alignment_tsv.py" "$out" "$readhashc_file" "$chunk_name"
        
        rm "$readhashc_file"

    else
        echo "Warning: No matching files found for pattern ${out}/${fastq}_chunk*readHashc.pkl"
        break  # Exit the loop if no files are found
    fi
done


# Combine the results from all chunks
python "${SCRIPT_DIR}/combine_BP.py" "$out" "$fastq"
python "${SCRIPT_DIR}/combine_readHashcc.py" "$out" "$fastq"
python "${SCRIPT_DIR}/combine_geneHashcc.py" "$out" "$fastq"
python "${SCRIPT_DIR}/combine_alignment_tsv.py" "$out" "$fastq"

# Generate breakpoint level output
python "${SCRIPT_DIR}/generate_BP.result.py" "$fastq" "$out" "$chimerabp"

# Generate transcript level output
python "${SCRIPT_DIR}/generate_ENSTresult.py" "$fastq" "$out"

# Generate gene level output
python "${SCRIPT_DIR}/generate_ENSGresult.py" "$fastq" "$out"

echo "Finished fusionBlaster"

if [ $? -eq 0 ]; then
    echo "Script executed successfully, removing files..."
    gzip -9 ${out}/${fastq}.readHashcc.pkl
    gzip -9 ${out}/${fastq}.readAlignment.tsv
else
    echo "Script failed or ran out of memory, files will not be removed."
fi