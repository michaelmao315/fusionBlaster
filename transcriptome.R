oldw <- getOption("warn")
options(warn = -1)
suppressMessages(library(ensembldb, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(plyr, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(AnnotationHub, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(GenomicFeatures, quietly = TRUE, warn.conflicts = FALSE))

genome_version <- commandArgs()[6]

if(genome_version == "grch37" || genome_version == "hg19"){ 
  suppressMessages(library(EnsDb.Hsapiens.v75, quietly = TRUE, warn.conflicts = FALSE))
  ensembledb <- EnsDb.Hsapiens.v75
  
  # Specify the path to the downloaded hg19.2bit file
  twoBitFile <- "/project/hlilab/software/fusionBlaster/hg19.2bit"
  dna <- rtracklayer::TwoBitFile(twoBitFile)
  
}else if(genome_version == "grch38" || genome_version == "hg38"){
  suppressMessages(library(EnsDb.Hsapiens.v86, quietly = TRUE, warn.conflicts = FALSE))
  ensembledb <- EnsDb.Hsapiens.v86
  dna <- ensembldb:::getGenomeTwoBitFile(ensembledb)
  
}else{
  print("invalid genome")
  quit(save="no")
}

# Processing ENSTs and extracting sequences
ensts <- suppressMessages(read.table(commandArgs()[7], header = FALSE, sep = '\t', col.names = c("ENST")))

ensts <- as.matrix(ensts)
for(i in 1:nrow(ensts)) {
  ENST <- ensts[i,1]
  Tx <- exonsBy(ensembledb, filter = TxIdFilter(ENST))
  
  # Modify the GRanges object (Tx) to add the "chr" prefix for GRCh37
  if(genome_version == "grch37" || genome_version == "hg19"){
    seqlevels(Tx) <- paste0("chr", seqlevels(Tx))
  }
  
  # Extract the transcript sequences using the modified chromosome names
  TxSeqs <- extractTranscriptSeqs(dna, Tx)
  vartmp <- as.data.frame(TxSeqs)
  name <- row.names(vartmp)
  Result <- paste(name, vartmp$x)
  cat(Result, sep = '\n')
}

options(warn = oldw)
