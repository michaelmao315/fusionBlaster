## position 1/commandArgs()[6] is the reference genome and must be grch37/hg19 or grch38/hg38
## position 2/commandArgs()[7] is the genomecoord.input file: must be a 4 column file with first column as up ENSG (ENSG00000169926.5), second column as down ENSG (ENSG00000169032.5), third column as 5'genebreakpoint (15:31658757:+), and fourth column as 3'genebreakpoint (15:66782056:+)
oldw <- getOption("warn")
options(warn = -1)
suppressMessages(library(ensembldb, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(plyr, quietly = TRUE, warn.conflicts = FALSE))

if(commandArgs()[6] == "grch37" || commandArgs()[6] == "hg19"){ 
  suppressMessages(library(EnsDb.Hsapiens.v75, quietly = TRUE, warn.conflicts = FALSE))
  ensembledb <- EnsDb.Hsapiens.v75
}else if(commandArgs()[6] == "grch38" || commandArgs()[6] == "hg38"){
  suppressMessages(library(EnsDb.Hsapiens.v86, quietly = TRUE, warn.conflicts = FALSE))
  ensembledb <- EnsDb.Hsapiens.v86
}else if(commandArgs()[6] == "grcm38" || commandArgs()[6] == "mm10"){
  suppressMessages(library(EnsDb.Mmusculus.v79, quietly = TRUE, warn.conflicts = FALSE))
  ensembledb <- EnsDb.Mmusculus.v79
}else{
  print("invalid genome")
  quit(save="no")
}

## create input file for script to read
input <- suppressMessages(read.table(commandArgs()[7], header = FALSE, sep = '\t', col.names = c("UpENSG", "DownENSG", "Up", "Down")))

## establish for loop                                                                                                                                     
for(i in 1:nrow(input)){                                                                                                                                  
upid <- input[i,1]
downid <- input[i,2]

## genome to Transcript for 5' and 3' parent                                                                                                              
upinput <- GRanges(input[i,3])
edbx <- filter(ensembledb, filter = ~ gene_id == upid)
upgnm <- genomeToTranscript(upinput, edbx)                                                                                                          
## if filtered database doesn't match anything, run unfiltered database
if(all(upgnm@unlistData@start == "-1")){
  upgnm <- genomeToTranscript(upinput, ensembledb)
}

downinput <- GRanges(input[i,4])
edbx <- filter(ensembledb, filter = ~ gene_id == downid)
downgnm <- genomeToTranscript(downinput, edbx)                                                                                                      
## if filtered database doesn't match anything, run unfiltered database
if(all(downgnm@unlistData@start == "-1")){
  downgnm <- genomeToTranscript(downinput, ensembledb)
}

## grab transcript coordinate with ENST, paste them all together, and filter out LRG entries
if(all(upgnm@unlistData@start != "-1")){
for(id in 1:length(upgnm@unlistData@start)){                                                                                                            
    if(grepl("LRG", upgnm@unlistData@NAMES[id]) == "FALSE" && exists("uptcoord") == "FALSE"){                                                         
      uptcoord <- paste(upgnm@unlistData@start[id], upgnm@unlistData@NAMES[id], sep = ",")                                                      
    }else if(grepl("LRG", upgnm@unlistData@NAMES[id]) == "FALSE" && exists("uptcoord") == "TRUE"){                                                    
      tmp <- paste(upgnm@unlistData@start[id], upgnm@unlistData@NAMES[id], sep = ",")                                                           
      uptcoord <- paste(uptcoord, tmp, sep = ":")                                                                                               
    }                                                                                                                                                 
} }else{ uptcoord <- "Intron Region"}
  
if(all(downgnm@unlistData@start != "-1")){
  for(id in 1:length(downgnm@unlistData@start)){                                                                                                        
    if(grepl("LRG", downgnm@unlistData@NAMES[id]) == "FALSE" && exists("downtcoord") == "FALSE"){                                                     
      downtcoord <- paste(downgnm@unlistData@start[id], downgnm@unlistData@NAMES[id], sep = ",")                                                
    }else if(grepl("LRG", downgnm@unlistData@NAMES[id]) == "FALSE" && exists("downtcoord") == "TRUE"){                                                
      tmp <- paste(downgnm@unlistData@start[id], downgnm@unlistData@NAMES[id], sep = ",")                                                       
      downtcoord <- paste(downtcoord, tmp, sep = ":")                                                                                           
    }                                                                                                                                                 
} }else { downtcoord <- "Intron Region"}

## create output file and bind all subsequent outputs                                                                                                     
if(i == 1){                                                                                                                                               
    Result <- data.frame(matrix(data = c(upid,downid,uptcoord,downtcoord), nrow = 1))                                                                              
}else{                                                                                                                                                  
  df <- data.frame(matrix(data = c(upid,downid,uptcoord,downtcoord), nrow = 1))                                                                                  
  Result <- rbind.fill(Result, df)                                                                                                                        
}                                                                                                                                                       

## remove variables to reset for LRG filter
rm(uptcoord)                                                                                                                                              
rm(downtcoord)

## close loop
}                                                                                                                                                         
rm(df)                                                                                                                                                    

write.table(Result, file = "", quote = FALSE, sep = "\t",                                                                                                 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,                                                                                          
            col.names = FALSE)                                                                                                                            

options(warn = oldw) 
