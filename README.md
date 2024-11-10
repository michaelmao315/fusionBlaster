Overview

FusionBlaster is a toolkit for generating a chimeric RNA reference database and quantifying chimeras from RNAseq data. This guide provides detailed instructions on installation, setup, and running analyses.

Installation

1. Clone the Repository: Clone the repository to get access to the necessary scripts:
git clone https://github.com/michaelmao315/fusionBlaster.git

After cloning, a folder named fusionBlaster containing all the necessary programs should be created.


2. Install Dependencies: Run the install.sh script to install all required dependencies:
./fusionBlaster/install.sh

Usage

Step 1: Create Reference Database
To create a reference database of chimeric RNAs, use the db.writer.sh script:
./fusionBlaster/db.writer.sh <chimeric_input_file> <genome_id> [output_directory]

Arguments:

<chimeric_input_file>: Path to the input file containing chimeric RNA information.

<genome_id>: Genome ID for the reference, currently supporting grch38 and grch37.

[output_directory] (optional): Output directory. By default, the script creates a directory in fusionBlaster/database to store the output files required for the quantitation program.

Input File Format: The input file should be a five-column file with the following format:  
<chimera_id>   <geneA_id>   <geneB_id>   <geneA_coordinates>   <geneB_coordinates>

Example: 

AAK1--RAB3IP   ENSG00000115977   ENSG00000127328   2:69460853:-   12:69794438:+

AARS--APBB2    ENSG00000090861   ENSG00000163697   16:70289421:-   4:40811474:+

Where:

Column 1: Chimera identifier.
Column 2: Ensemble gene ID for gene A.
Column 3: Ensemble gene ID for gene B.
Column 4: Genomic coordinates of gene A (format: chromosome:position
).
Column 5: Genomic coordinates of gene B (format: chromosome:position
).

Step 2: Quantify Chimeras

Once the reference database is created, use the fusionBlaster.sh script to quantify each chimera from RNAseq data:

./fusionBlaster/fusionBlaster.sh <data_dir> <output_dir> <SRR_id> <refDB.fa> <filtered_output.tsv> <num_cores>

Arguments:

<data_dir>: Directory containing paired-end sequencing data files (e.g., SRR8615300_1.fq.gz, SRR8615300_2.fq.gz).

<output_dir>: Directory to store the output.

<SRR_id>: Identifier for RNAseq files, used to locate paired-end files.

<refDB.fa>: Reference database file generated in the previous step.

<filtered_output.tsv>: Filtered output file containing transcript information for each chimera, also generated in the previous step, should be in the same directory as the reference databse.

<num_cores>: Number of cores to use for parallel processing.

Example command

./fusionBlaster/fusionBlaster.sh ./data/SRR8615300 ./output/SRR8615300_3_output SRR8615300 ./fusionBlaster/database/refDB.fa ./fusionBlaster/database/filtered_output.tsv 2


In this example:

./data/SRR8615300: Directory with sequencing data (files named SRR8615300_1.fq.gz and SRR8615300_2.fq.gz).
./output/SRR8615300_3_output: Output directory.
SRR8615300: SRR ID used to locate RNAseq files.
./fusionBlaster/database/refDB.fa: Reference database file.
./fusionBlaster/database/filtered_output.tsv: Filtered chimeric input file.
2: Number of cores for parallel processing.


