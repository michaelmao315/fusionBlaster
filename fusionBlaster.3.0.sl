#!/bin/bash
#SBATCH --account=hlilab
#SBATCH --partition=standard
#SBATCH --time=0-95:59:59
#SBATCH --nodes=1

in=$1
out=$2
fastq=$3
refDB=$4
chimerabp=$5
cores=$6

print_elapsed_time() {
    end_time=$(date +%s)
    elapsed_time=$((end_time - start_time))
    echo "Elapsed time: $(date -ud "@$elapsed_time" +'%H:%M:%S')"
}

start_time=$(date +%s)

if [[ $in == "" || $out == "" || $fastq == "" || refDB == "" ]]; then 
echo "Proper Usage:"
echo 'sbatch --error=[uniq].err --output=[uniq].out --job-name=[uniq] fusionBlaster.sl /full/path/to/directory/of/fastqs /full/path/to/output/directory uniqFASTQid full/path/to/refDB.fa /full/path/to/chimeras.input.tsv'
echo '=================================================================='
echo '[uniq] needs to be unique to each set of fastq files; should match uniqFASTQid'
echo '/full/path/to/directory/of/fastqs should be full path to a central directory where all fastq files to be analyzed are located'
echo '/full/path/to/output/directory should be single output directory where results from all fastqs anlayzed by same chimeras.input.tsv are located'
echo 'uniqFASTQid should be used from fastq files as [uniqFASTQid]_R1.fq.gz or [uniqFASTQid]_2.fastq.gz'
echo '/full/path/to/refDB.fa should be full path to the refDB.fa file created as STDOUT by the refDBwriter.sl'
echo '/full/path/to/chimeras.input.tsv is the full pth to the original input file with chimera names, ENSGs, and breakpoint info'
exit
fi

module load anaconda

fq1=$(ls $in | grep "^$fastq"_ | grep -f <( echo fq ; echo fastq) |  sort -u | head -1)
fq2=$(ls $in | grep "^$fastq"_ | grep -f <( echo fq ; echo fastq) |  sort -u | tail -1)
if [[ $fq1 == "" || $fq2 == "" || $fq1 != *".gz" || $fq2 != *".gz" ]]; then 
echo 'fastq files were either not found in the specified input directory or are not named correctly'
echo 'fastq files must be paired end fastq files name like:'
echo 'uniq_1.fastq.gz uniq_2.fastq.gz or uniq_1.fq.gz uniq_2.fq.gz'
echo 'uniq_R1.fastq.gz uniq_R2.fastq.gz or uniq_R1.fq.gz uniq_R2.fq.gz'
exit
fi

echo ${in}/$fq1 ${in}/$fq2

readLength=$(for read in $(zcat ${in}/$fq1 | tail -n +2 | awk 'NR%4==1' | head -2000); do echo $read | wc -m ; done | awk '{ sum += $1-1 } END { if (NR > 0) print sum / NR }' | cut -d. -f1)
cdir=$(pwd)

echo $readLength
echo $cdir
echo $out

source ~/.bashrc
if [ ! -d "$out" ]; then
  mkdir -p "$out"
fi

#convert fastqs to a fasta with read_1 and read_2 pasted together
zcat ${in}/$fq1 | paste - - - - | cut -f1,2 | tr "@""\t" ">""\n" > ${out}/${fastq}_tmp1
paste <(zcat ${in}/$fq2 | paste - - - - | cut -f1 | tr "@" ">") <(zcat ${in}/$fq2 | paste - - - - | cut -f2 | tr ACTG TGAC | rev) | tr "\t" "\n" > ${out}/${fastq}_tmp2
paste ${out}/${fastq}_tmp1 ${out}/${fastq}_tmp2 | awk -v IFS="\t" -v OFS="" '{if ($1~">") print $1; else print $1,$2 }' > ${out}/${fastq}_query.fa
rm -rf ${out}/${fastq}_tmp1 ${out}/${fastq}_tmp2

echo "Finished everything before pBlat"
print_elapsed_time


start_time=$(date +%s)

#run pBLAT and filter alignments

split -l 1000000  ${out}/${fastq}_query.fa ${out}/${fastq}_chunk_
rm ${out}/${fastq}_query.fa

# Main processing loop
cd $out

for chunk in ${fastq}_chunk_*; do
    chunk_prefix=$(basename "$chunk" .fa)
    
    echo "Processing chunk: $chunk_prefix"

    # Run pBLAT on the chunk
    conda activate pblat
    pblat -fastMap -threads=$cores -out=blast8 $refDB $chunk ${chunk_prefix}_pblat.out
    conda deactivate

    # Filter
    awk -v OFS='\t' -v L=$readLength '{ if ($3>98 && $4>L-3) print $1,$2,$7,$8,$9,$10 }' ${chunk_prefix}_pblat.out > ${chunk_prefix}_blast.tsv
    
    rm $chunk ${chunk_prefix}_pblat.out
    
    # Ensure every chunk of blast.tsv is below 10000000*cores lines
    lines=$(wc -l < ${chunk_prefix}_blast.tsv)
    if [ $lines -gt $((10000000 * cores)) ]; then
        python ./fusionBlaster/split_file.py ${chunk_prefix}_blast.tsv ${chunk_prefix}_blast_division_ $((10000000 * cores))
        rm ${chunk_prefix}_blast.tsv
    else
        mv ${chunk_prefix}_blast.tsv ${chunk_prefix}_blast_division_0001.tsv
    fi

    # Run fusionBlaster_1.0.py on all blast TSV divisions
    for file in ${chunk_prefix}_blast_division_*.tsv; do
        division=$(basename "$file" .tsv | rev | cut -d '_' -f 1 | rev)
        python ./fusionBlaster/fusionBlaster.2.0.py $file ${chunk_prefix}_${division} $out
    done
    
    # Remove temporary files
    rm ${chunk_prefix}_blast_division_*.tsv
    
    echo "Finished processing chunk: $chunk_prefix"
done

echo "All chunks processed"

echo "Finished pBlat Loop"
print_elapsed_time


start_time=$(date +%s)

for readhashc_file in ${out}/${fastq}_chunk*readHashc.pkl; do
    if [ -f "$readhashc_file" ]; then
        chunk_name=$(basename "$readhashc_file" .readHashc.pkl)
        
        # Run count_BP.py for each chunk
        python ./fusionBlaster/count_BP.py "$out" "$chimerabp" "$readhashc_file" "$chunk_name"

        # Run count_readHashcc.py for each chunk
        python ./fusionBlaster/count_readHashcc.py "$out" "$readhashc_file" "$chunk_name"

        # Run count_geneHashcc.py for each chunk
        python ./fusionBlaster/count_geneHashcc.py "$out" "$readhashc_file" "$chunk_name"

        # Generate alignment TSV for each chunk
        python ./fusionBlaster/generate_alignment_tsv.py "$out" "$readhashc_file" "$chunk_name"
        
        rm "$readhashc_file"

    else
        echo "Warning: No matching files found for pattern ${out}/${fastq}_chunk*readHashc.pkl"
        break  # Exit the loop if no files are found
    fi
done

start_time=$(date +%s)
cd $out

# Combine the results from all chunks
python ./fusionBlaster/combine_BP.py $out $fastq
python ./fusionBlaster/combine_geneHashcc.py $out $fastq
python ./fusionBlaster/combine_alignment_tsv.py "$out" "$fastq"

# Generate breakpoint level output
python ./fusionBlaster/generate_BP.result.py ${fastq} $out $chimerabp 

# Generate transcript level output
python ./fusionBlaster/generate_ENSTresult.py ${fastq} $out

# Generate gene level output
python ./fusionBlaster/generate_ENSGresult.py ${fastq} $out


echo "Finished fusionBlaster"
print_elapsed_time

if [ $? -eq 0 ]; then
    echo "Script executed successfully, removing files..."
    gzip -9 ${out}/${fastq}.readHashcc.pkl
    gzip -9 ${out}/${fastq}.readAlignment.tsv
else
    echo "Script failed or ran out of memory, files will not be removed."
fi 