library(magrittr)
library(biomaRt)
library(R6)

# Placeholder for maybe_create_directory
maybe_create_directory <- function(filepath) {
  if (!dir.exists(filepath)) {
    dir.create(filepath, recursive = TRUE)
  }
}

# Placeholder for log_message
log_message <- function(message) {
  cat(paste0("[INFO] ", message, "\n"))
}

#' BioMartData R6 class
#'
#' An R6 class to fetch, manipulate, and save gene data from the Ensembl BioMart 
#' database. It allows fetching data for specific chromosomes and datasets, 
#' and stores all the fetched data internally for further processing.
#'
#' @field ensembl The Ensembl BioMart object used for fetching data.
#' @field dataset The dataset from which data is being fetched.
#' @field chromosomes The list of chromosomes for which data will be fetched.
#' @field combined_data A combined data frame of all the fetched data.
#'
#' @description 
#' The `initialize` method sets up the Ensembl BioMart object based on the provided `biomart` 
#' and `dataset` parameters. It also initializes the `chromosomes` field with the list of 
#' chromosomes available in the dataset.
#'
#' The `get_data` method fetches gene data for specified chromosomes or for all available 
#' chromosomes if none are specified. It also manages directory creation, error handling, 
#' data combination, and temporary file handling.
#'
#' @examples
#' # Create an instance of BioMartData
#' biomart <- BioMartData$new(biomart = "ensembl", dataset = "btaurus_gene_ensembl")
#' 
#' # Retrieve and process the data
#' biomart$get_data(chromosomes = c("1", "2"), filepath = ".")
#' 
#' @export
BioMartData <- R6Class(
  "BioMartData",
  public = list(
    ensembl = NULL,
    dataset = NULL,
    chromosomes = NULL,
    combined_data = NULL,

    #' @description Initialize BioMartData object.
    #' @param biomart Character string. The biomart to use.
    #' @param dataset Character string. The dataset to use.
    initialize = function(biomart, dataset) {
      self$ensembl <- useEnsembl(biomart = biomart)
      self$ensembl <- useDataset(dataset = dataset, mart = self$ensembl)
      self$dataset <- dataset
      self$chromosomes <- getBM(attributes = 'chromosome_name', mart = self$ensembl)[, 1]
    },

    #' @description Fetch gene data from Ensembl BioMart.
    #' @param chromosomes Character vector (default NULL). The chromosomes for which data will be fetched.
    #' @param filepath Character string (default PROTEOMES_PATH). Path where the fetched data will be stored.
    get_data = function(chromosomes = NULL, filepath = PROTEOMES_PATH) {
      maybe_create_directory(filepath)
      
      if (!is.null(chromosomes)) {
        chromosomes <- intersect(chromosomes, self$chromosomes)
      } else {
        chromosomes <- self$chromosomes
      }
      
      log_message("Fetching data for chromosomes...")
      data_list <- list()
      file_list <- list()  # Create a list to hold the filenames
      
      for (chrom in chromosomes) {
        log_message(paste("Fetching data for chromosome", chrom))
        tryCatch({
          annotLookup <- getBM(
            filters = c("chromosome_name"),
            values = list(chrom),
            mart = self$ensembl,
            attributes = c(
              'external_gene_name',
              'uniprot_gn_id',
              'chromosome_name',
              'start_position',
              'end_position',
              'peptide'
            ),
            uniqueRows = TRUE
          )
          
          annotLookup <- subset(annotLookup, uniprot_gn_id != '')
          filename <- file.path(filepath, paste0("annotLookup_", self$dataset, "_", chrom, ".csv"))
          write.csv(annotLookup, filename, row.names = FALSE)
          data_list[[chrom]] <- annotLookup
          file_list <- c(file_list, filename)  # Add the filename to the list
        },
        error = function(e) {
          log_message(paste("Failed to fetch data for chromosome", chrom, "due to error:", e))
        })
      }
      
      log_message("Data fetching completed.")
      
      self$combined_data <- do.call(rbind, data_list)
      self$combined_data <- self$combined_data[!(self$combined_data$peptide == "Sequence unavailable"), ]
      
      self$combined_data <- self$combined_data[, c("external_gene_name", "uniprot_gn_id", "chromosome_name", "start_position", "end_position", "peptide")]
      
      colnames(self$combined_data) <- c("GeneName", "SwissID", "chromosome", "start", "end", "seq")
      
      combined_filename <- file.path(filepath, paste0(self$dataset, "_combined.csv"))
      file_list <- c(file_list, combined_filename)
      write.csv(self$combined_data, combined_filename, row.names = FALSE)
      
      self$combined_data <- lapply(self$combined_data, function(column) {
        quoted_values <- paste0('"', column, '"')
        return(quoted_values)
      })
      
      quoted_filename <- file.path(filepath, paste0(self$dataset, ".csv"))
      write.csv(self$combined_data, quoted_filename, row.names = FALSE, quote = FALSE)
      
      # Remove the temporary files
      log_message(paste("Removing temporary files:", file_list))
      
      full_file_paths <- normalizePath(file.path(file_list), mustWork = FALSE)
      file.remove(full_file_paths)
      
      log_message("Downloading proteome completed.")
    }
  )
)

# Create an instance of BioMartData
# biomart <- BioMartData$new(biomart = "ensembl", dataset = "btaurus_gene_ensembl")

# Retrieve and process the data
# biomart$get_data()

# To download specific chromosomes, use the following:
# biomart$get_data(chromosomes = c("1", "2"))
