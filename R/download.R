require('biomaRt')

BioMartData <- R6::R6Class(
  "BioMartData",
  public = list(
    ensembl = NULL,
    dataset = NULL,
    chromosomes = NULL,
    initialize = function(dataset, chromosomes) {
      self$ensembl <- useEnsembl(biomart = "genes")
      self$dataset <- useDataset(dataset = dataset, mart = self$ensembl)
      self$chromosomes <- chromosomes
    },
    get_data = function(filepath) {
      for(chrom in self$chromosomes) {
        annotLookup <- getBM(
          filters=c("chromosome_name"),
          values=list(chrom),
          mart = self$dataset,
          attributes = c(
            'external_gene_name',
            'uniprot_gn_id',
            'chromosome_name',
            'start_position',
            'end_position',
            'peptide'),
          uniqueRows=TRUE)
        
        annotLookup <- subset(annotLookup, uniprot_gn_id != '')
        write.csv(annotLookup, file.path(filepath, paste0("annotLookup_", chrom, ".csv")), row.names = FALSE)
      }
    },
    combine_csv = function(filepath) {
      files <- list.files(filepath, pattern = ".csv$")
      data_list <- lapply(files, function(file) read.csv(file.path(filepath, file)))
      combined_data <- do.call(rbind, data_list)
      combined_data <- combined_data[!(combined_data$peptide=="Sequence unavailable"),]
      write.csv(combined_data, file.path(filepath, "combined.csv"), row.names = FALSE)
      return(combined_data)
    },
    add_quotes = function(data, output_file) {
      quoted_data <- lapply(data, function(column) {
        quoted_values <- paste0('"', column, '"')
        return(quoted_values)
      })
      write.csv(quoted_data, output_file, row.names = FALSE, quote = FALSE)
    }
  )
)

# Usage:
# chromosomes <- c(as.character(1:29), "X")
# data <- BioMartData$new(dataset = "btaurus_gene_ensembl", chromosomes = chromosomes)
# data$get_data(filepath = "C:/Projectfiletest/btaurus_proteome_Uniprot")
# combined <- data$combine_csv("C:/Projectfiletest/btaurus_proteome_Uniprot")
# data$add_quotes(combined, "C:/Projectfiletest/btaurus_proteome_Uniprot/combined_quoted.csv")
