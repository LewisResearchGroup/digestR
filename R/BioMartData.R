library(magrittr)

BioMartData <- R6::R6Class(
  "BioMartData",
  public = list(
    ensembl = NULL,
    dataset = NULL,
    chromosomes = NULL,
    combined_data = NULL,
    initialize = function(biomart, dataset) {
      self$ensembl <- useEnsembl(biomart = biomart)
      self$ensembl <- useDataset(dataset = dataset, mart = self$ensembl)
      self$dataset <- dataset
      self$chromosomes <- getBM(attributes = 'chromosome_name', mart = self$ensembl)[,1]
    },
    get_data = function(chromosomes = NULL) {
      filepath <- "data"
      
      if (!is.null(chromosomes)) {
        chromosomes <- intersect(chromosomes, self$chromosomes)
      } else {
        chromosomes <- self$chromosomes
      }

      if (!dir.exists(filepath)) {
        dir.create(filepath)
      }

      print(getwd())
      print("Fetching data for chromosomes...")
      data_list <- list()
      file_list <- list()  # Create a list to hold the filenames
      
      for (chrom in chromosomes) {
        print(paste("Fetching data for chromosome", chrom))
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
              'peptide'),
            uniqueRows = TRUE
          )

          annotLookup <- subset(annotLookup, uniprot_gn_id != '')
          filename <- file.path(filepath, paste0("annotLookup_", self$dataset, "_", chrom, ".csv"))
          write.csv(annotLookup, filename, row.names = FALSE)
          data_list[[chrom]] <- annotLookup
          file_list <- c(file_list, filename)  # Add the filename to the list
        },
        error = function(e) {
          print(paste("Failed to fetch data for chromosome", chrom, "due to error:", e))
        })
      }
      
      print("Data fetching completed.")
      
      self$combined_data <- do.call(rbind, data_list)
      self$combined_data <- self$combined_data[!(self$combined_data$peptide == "Sequence unavailable"), ]
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
      print(paste("Removing temporary files:", file_list))

      full_file_paths <- normalizePath(file.path(file_list), mustWork = FALSE)
      file.remove(full_file_paths)
      
      print("Data processing completed.")
    }

  )
)

# Create an instance of BioMartData
# data <- BioMartData$new(biomart = "ensembl", dataset = "btaurus_gene_ensembl")

# Retrieve and process the data
# data$get_data()

# To download specific chromosomes, use the following:
# data$get_data(chromosomes = c("1", "2"))
