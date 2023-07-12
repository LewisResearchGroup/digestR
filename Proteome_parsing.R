

##################################################################################

## Parsing uniprot fasta files
## order of data in file:   GeneName, UniProtKB, Chromosome, GeneStart, GeneEnd, Peptide
bostaurus <- paste(readLines( myOpen(".txt", list(txt = "text file"), multiple = FALSE), warn = FALSE), collapse="")


#matches <- cpp_gregexpr2(targetSeq, pep$sequence[i])
genes <- cpp_gregexpr2(bostaurus, '>')
seqEndIndices <- (genes[2:length(genes)] - 1)
seqEndIndices <- c(seqEndIndices, nchar(bostaurus))
seps <- cpp_gregexpr2(bostaurus, '|')

btdata <- data.frame(geneName = character(), uniProtKB = character(), chromosome = character(), 
                     start=character(), end=character(), seq=character(), stringsAsFactors=FALSE)

for ( i in 1:length(genes)){
  
  btdataEntry <- data.frame(geneName = '', uniProtKB = '', chromosome = '', 
                            start='', end='', seq='', stringsAsFactors=FALSE)
  
  sLine <- substr(bostaurus, genes[i], seqEndIndices[i])
  
  idx <- regexpr('GN=', sLine)[[1]][1]
  if(idx != -1)
  {
    idx <- idx + 3
    idx2 <- regexpr(' ', substr(sLine, idx, nchar(sLine)))[[1]][1]
    
    if(idx2 != -1)
    {
      idx2 <- idx + idx2 - 2 # subtract 1 to get rid of space, subtract another to account to make index right
      btdataEntry$geneName <- substr(sLine, idx, idx2)
    }
    else
      btdataEntry$geneName <- ''	
  }else
    btdataEntry$geneName <- ''
  
  seps <- cpp_gregexpr2(sLine, '|')
  
  if(length(seps)>1 & (seps[1] < seps[2]))
    btdataEntry$uniProtKB <- substr(sLine, seps[1] + 1, seps[2] - 1)
  else
    btdataEntry$uniProtKB <- '' 
  
  idx <- regexpr('SV=', sLine)[[1]][1]
  if(idx != -1)
  {
    idx <- idx + 3
    idx2 <- regexpr('[a-zA-Z]', substr(sLine, idx, nchar(sLine)))[[1]][1]
    
    if(idx2 != -1)
    {
      idx2 <- idx + idx2 - 1 # subtract 1 to get rid of space, subtract another to account to make index right
      btdataEntry$seq <- substr(sLine, idx2, nchar(sLine))
    }
    else
      btdataEntry$seq <- ''	
  }else
    btdataEntry$seq <- ''
  
  if(nchar(btdataEntry$seq)>0)
    btdata <- rbind(btdata, btdataEntry)
  
}


##################################################################################

## Parsing biomart fasta files

bostaurus <- paste(readLines( myOpen(".txt", list(txt = "text file"), multiple = FALSE), warn = FALSE), collapse="")

bostaurus <- gsub("_coding", "_coding:", bostaurus)
bostaurus <- gsub("transcript", ":transcript", bostaurus)

genes <- cpp_gregexpr2(bostaurus, '>')               # split your txt file into string that start at > and end at >

# create a map of where the > character is in your file we start at 2 because we know where the first index start 
seqEndIndices <- (genes[2:length(genes)] - 1)   
# append to your data the valume of nchar() -> nchar = number of characters in bostarus
seqEndIndices <- c(seqEndIndices, nchar(bostaurus)) 
offset <- 4

indices <- vector(mode = "list", length = 5)
indices <- 1:5

names(indices) <- c('uniProtKB', 'chrom', 'geneStart', 'geneEnd','seq')

btdata <- data.frame(geneName = character(), uniProtKB = character(), chromosome = character(), 
                      start=character(), end=character(), seq=character(),stringsAsFactors=FALSE)

btdata_final <- data.frame(geneName = '', uniProtKB = '', chromosome = '', 
                          start='', end='', seq='', stringsAsFactors=FALSE)
# length(genes)
for ( i in 1:length(genes)) {
  
  btdataEntry <- data.frame(geneName = '', uniProtKB = '', chromosome = '', 
                            start='', end='', seq='', stringsAsFactors=FALSE)
  
  sLine <- substr(bostaurus, genes[i], seqEndIndices[i])
  
  seps <- cpp_gregexpr2(sLine, ':')                 # returns a list of positions where your character ':' exist in your string

   if((seps[indices['uniProtKB']]-1)>=2)
     btdataEntry$geneName <- substr(sLine, 2, seps[indices['uniProtKB']]-1)
   else
     btdataEntry$geneName <- ''
  
  if(seps[indices['uniProtKB']]+1 <= seps[indices['chrom']]-1){
     btdataEntry$uniProtKB <- substr(sLine, seps[indices['uniProtKB']]+1, seps[indices['chrom']]-1) # modfied to store in btdataentry
   }else
     btdataEntry$uniProtKB <- '' 
  
   if(seps[indices['chrom']]+1 <= seps[indices['geneStart']]-1)
     btdataEntry$chromosome <- substr(sLine, seps[indices['chrom']]+1,  seps[indices['geneStart']]-1)
   else
     btdataEntry$chromosome <- 'ukn'
  
   if(seps[indices['geneStart']]+1 <= seps[indices['geneEnd']]-1)
     btdataEntry$start <- substr(sLine, seps[indices['geneStart']]+1,  seps[indices['geneEnd']]-1)
   else
     btdataEntry$start <- ''
  
  if(seps[indices['geneEnd']]+1 <= seps[indices['seq']]-1)
     btdataEntry$end <- substr(sLine, seps[indices['geneEnd']]+1,  seps[indices['seq']]-1)
   else
     btdataEntry$end <- ''
  
  geneSeq <- strsplit(sLine,'coding:')
  btdataEntry$seq<- geneSeq[[1]][2]
    
    btdata_final <- rbind(btdata_final,btdataEntry)
  #Could change the data set with strsplit function strsplit(sLine, ":")
}
write.csv(btdata_final,file.choose(),row.names = FALSE) 

##################################################################################

## Parsing biomart fasta files with Ensembl IDs

bostaurus <- paste(readLines( myOpen(".txt", list(txt = "text file"), multiple = FALSE), warn = FALSE), collapse="")

bostaurus <- gsub("_gene", "_gene:", bostaurus)
bostaurus <- gsub("_coding", "_coding:", bostaurus)
bostaurus <- gsub("]", "]:", bostaurus)
bostaurus <- gsub(" description", ":description", bostaurus)

genes <- cpp_gregexpr2(bostaurus, '>')
seps <- unlist(strsplit(bostaurus,':'))

ensmblID <- seps[[1]][1]
geneName <- seps[1][13]
chrom <- seps[[1]][3]
geneStart <- seps[[1]][4]
geneEnd <- seps[[1]][5]

ensmblID <- genes[which(genes == ">")+1]
geneName <- genes[which(genes == ">")+1]
proE <- ogenes[which(genes == ">")+1]
pro <- out[which(out == "SV=")+1]

btdata <- data.frame(ensmblID = character(), geneName = character(), chromosome = character(), 
                     start=character(), end=character(), seq=character(),stringsAsFactors=FALSE)


################################################################################################
## Parsing biomart fasta files with Ensembl IDs

bostaurus <- paste(readLines( myOpen(".txt", list(txt = "text file"), multiple = FALSE), warn = FALSE), collapse="")

bostaurus <- gsub("_gene", "_gene:", bostaurus)
bostaurus <- gsub("_coding", "_coding:", bostaurus)
bostaurus <- gsub("]", "]:", bostaurus)
bostaurus <- gsub(" description", ":description", bostaurus)

genes <- cpp_gregexpr2(bostaurus, '>')

seqEndIndices <- (genes[2:length(genes)] - 1)   

seqEndIndices <- c(seqEndIndices, nchar(bostaurus)) 

indices <- vector(mode = "list", length = 5)
indices <- 1:5

names(indices) <- c('uniProtKB', 'chrom', 'geneStart', 'geneEnd','seq')

btdata <- data.frame(ensemblID = character(), geneName = character (), uniProtKB = character(), chromosome = character(), 
                     start=character(), end=character (), seq=character(),stringsAsFactors=FALSE)

btdata_final <- data.frame(ensemblID = '', geneName = '', uniprotKB = '', chromosome = '', 
                           start='', end='', seq='', stringsAsFactors=FALSE)
#length(genes)
for ( i in 1:length(genes)){
  
  btdataEntry <- data.frame(ensemblID = '', geneName = '', uniProtKB = '', chromosome = '', 
                            start='', end='', seq='', stringsAsFactors=FALSE)
  
  sLine <- substr(bostaurus, genes[i], seqEndIndices[i])
  
  seps <- cpp_gregexpr2(sLine, ':') 
  
  
  if((seps[indices['uniProtKB']]-1)>=2)
    btdataEntry$ensemblID <- substr(sLine, 2, seps[indices['uniProtKB']]-24)
  else
    btdataEntry$ensemblID <- ''
  
  if(seps[indices['uniProtKB']]+1 <= seps[indices['chrom']]-1){
    btdataEntry$uniProtKB <- substr(sLine, seps[indices['uniProtKB']]+1, seps[indices['chrom']]-1)
  }else
    btdataEntry$uniProtKB <- '' 
  
  if(seps[indices['chrom']]+1 <= seps[indices['geneStart']]-1)
    btdataEntry$chromosome <- substr(sLine, seps[indices['chrom']]+1,  seps[indices['geneStart']]-1)
  else
    btdataEntry$chromosome <- 'ukn'
  
  if(seps[indices['geneStart']]+1 <= seps[indices['geneEnd']]-1)
    btdataEntry$start <- substr(sLine, seps[indices['geneStart']]+1,  seps[indices['geneEnd']]-1)
  else
    btdataEntry$start <- ''
  
  if(seps[indices['geneEnd']]+1 <= seps[indices['seq']]-1)
    btdataEntry$end <- substr(sLine, seps[indices['geneEnd']]+1,  seps[indices['seq']]-1)
  else
    btdataEntry$end <- ''
  
  splits <- strsplit(sLine,":")
  geneName <- splits[[1]][13]
  
  btdata_final <- rbind(btdata_final,btdataEntry)
  
}

# write.csv(btdata_final,file.choose(),row.names = FALSE) 

#########################################################################################################################

#Parsing and match using the biomaRt package
#Example using bostaurus 

require('biomaRt')

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
head(datasets)a
searchDatasets(mart = ensembl, pattern = "btaurus")
#21 btaurus_gene_ensembl Cow genes (ARS-UCD1.2) ARS-UCD1.2
ensembl <- useDataset(dataset = "btaurus_gene_ensembl", mart = ensembl)

#filters = listFilters(ensembl)
#attributes = listAttributes(ensembl)
  #Searching peptide sequence is too long for biomart and timeout
    #search by chromosomes

chrom = c(1,2)
annotLookup <- getBM(
  filters=c("chromosome_name"),
  values=list(chrom),
  mart = ensembl,
  attributes = c(
    'external_gene_name',
    'uniprot_gn_id',
    'chromosome_name',
    'start_position',
    'end_position',
    'peptide'),
  uniqueRows=TRUE)

head(subset(annotLookup, uniprot_gn_id != ''), 20)[,-4]
write.csv(annotLookup,file.choose(),row.names = FALSE)

# chrom = c(1,2,"X")
# annotLookup <- getBM(
#   filters=c("chromosome_name"),
#   values=list(chrom),
#   mart = ensembl,
# attributes = c(
#   'ensembl_gene_id',
#   'ensembl_peptide_id',
#   'uniprot_gn_id',
#   'chromosome_name',
#   'start_position',
#   'end_position',
#   'peptide'),
# uniqueRows=TRUE)


# Combine .csv files
filepath <- 'C:/Projectfiletest/btaurus_proteome_Uniprot'
setwd(filepath)

file_list <- list.files(filepath, pattern=NULL, all.files=FALSE,
                        full.names=FALSE)

file_path_full <- paste(filepath,'/',file_list[1],sep = "")
f <- read.csv(file =file_path_full, header = TRUE)

for (i in 2:length(file_list)){
  file_path_full <- paste(filepath,'/',file_list[i],sep = "")
  f<-rbind(f, read.csv(file =file_path_full, header = TRUE))
}
#write.csv(f,file.choose(),row.names = FALSE)

##Remove genes that do not possess known protein sequences 

fclean <- f[!(f$peptide=="Sequence unavailable"),]
write.csv(fclean,file.choose(),row.names = FALSE)
################################################################################
## Add quotation marks back in

# Read the original CSV file

data <- read.csv(file.choose(),header=T)

# Specify the delimiter used in the CSV file
delimiter <- ","

# Add quotation marks to each value
quoted_data <- lapply(data, function(column) {
  quoted_values <- paste0('"', column, '"')
  return(quoted_values)
})

write.csv(quoted_data, file.choose(), row.names = FALSE, quote = FALSE)

# Write the modified data to a new CSV file
output_file <- readline("Enter the path for the modified CSV file: ")
write.csv(quoted_data, output_file, row.names = FALSE, quote = FALSE)


# Alternatively, if you want to overwrite the original file, use the following:
# write.csv(quoted_data, original_file, row.names = FALSE, quote = FALSE)


#Change the headers name to:
#
################################################################################

# To upload the proteome on DigestR you need to:

# add the proteome to DigestR species list and species files in 2017_09_16_rNMR_TravisR
### Search for: speciesList and speciesFiles in the global environment
### i.e. speciesList = c('Homo sapiens', 'Plasmodium falciparum (3D7)', 'Pseudomonas aeruginosa 01', 
#                 'Pseudomonas aeruginosa 14',), speciesFiles = c('chrHs.csv', 'chr3D7.csv', 'pa01.csv', 'pa14.csv',),

#add loadSpecies in digestomics_MK4.R

# Modify and loadgene functions based on the following model: 
# initTaurus <- function(sFileName = '')
# {
#   if (sFileName == '')
#     sFileName <- paste(getwd(), '/chrBtaurus1.csv', collapse='', sep='')
#   
#   taurus <- loadTaurusGenes(sFileName)
#   
#   return(taurus)
# }
# 
# #	loadBosTaurusGenes()
# loadTaurusGenes <- function(sTaurusFilename = "")
# {
#   ## Read Human Proteome to Genome Map file
#   if (sTaurusFilename == "")
#   {
#     dfTaurus <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
#   }else{
#     dfTaurus <- read.csv( sTaurusFilename, head = TRUE, stringsAsFactors = FALSE)
#   }
#   
#   fileInfo <- file.info(sTaurusFilename)
#   TaurusID <- as.integer(fileInfo$mtime)
#   
#   ## rename "chromosome" to "chrom" for consistency and brevity
#   names(dfTaurus)[names(dfTaurus) == 'chromosome'] <- "chrom"
#   
#   ## make sure all x chromosomes are in same case
#   dfTaurus$chrom[which(dfTaurus$chrom == "X")] <- "x"
#   
#   ## make sure unknown designations share a tag
#   dfTaurus$GeneName[dfTaurus$GeneName == ""] <- "UKN"
#   dfTaurus$SwissID[dfTaurus$SwissID == ""] <- "UKN"
#   dfTaurus <- subset(dfTaurus, select = c("GeneName", "seq", "chrom", "start"))
#   names(dfTaurus)[1] <- "name"
#   
#   return(prepareSpecies("BosTaurus", TaurusID, dfTaurus))
# }


# if The species list do not update automatically so when loading DigestR load the species list in a separate command => See SpeciesListUpdate.R
