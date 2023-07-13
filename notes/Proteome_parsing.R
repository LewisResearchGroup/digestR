#Parsing proteomes using the biomaRt package
 
#To install BiomaRt, enter the following lines: 
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("biomaRt")

#Load biomaRt 
  require('biomaRt')

#BiomaRt example using btaurus

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

  chrom = c(1,2, "X")
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
write.csv(f,file.choose(),row.names = FALSE)

#Remove genes that do not possess known protein sequences 

fclean <- f[!(f$peptide=="Sequence unavailable"),]
write.csv(fclean,file.choose(),row.names = FALSE)

## Add quotation marks back in

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


# To upload the proteome on DigestR you need to:

# In 2017_09_16_rNMR_Travis_DD_Edits.R

# Search for: speciesList and speciesFiles in the global environment, name the species and add the .csv file
# i.e. speciesList = c('Homo sapiens', 'Plasmodium falciparum (3D7)', 'Pseudomonas aeruginosa 01', 
#                 'Pseudomonas aeruginosa 14',), speciesFiles = c('chrHs.csv', 'chr3D7.csv', 'pa01.csv', 'pa14.csv',),

# In digestR_Code_Dimitri_final.R
# add the csv file to the loadspecies function 

# Modify loadgene functions based on the following model: 
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
# #	loadTaurusGenes()
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
