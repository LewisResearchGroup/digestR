library(tcltk)

biomartDownloadWindow <- function() {
  # Create a new top-level window
  tt <- tktoplevel()
  tkwm.title(tt, "DigestR BioMart Downloader")
  tkwm.geometry(tt, "400x200")  # Set the window size

  # Tcl variables for entry widgets
  biomartVar <- tclVar("genes")
  datasetVar <- tclVar("btaurus_gene_ensembl")
  chromosomesVar <- tclVar("1, 2")

  # Create and position the controls
  lab1 <- tklabel(tt, text = "Biomart")
  ent1 <- tkentry(tt, textvariable = biomartVar, width=30)
  
  lab2 <- tklabel(tt, text = "Dataset")
  ent2 <- tkentry(tt, textvariable = datasetVar, width=30)
  
  lab3 <- tklabel(tt, text = "Chromosomes")
  ent3 <- tkentry(tt, textvariable = chromosomesVar, width=30)

  # When the button is clicked, this function will be executed
  onDownloadClick <- function() {
    mart <- tclvalue(biomartVar)
    dataset <- tclvalue(datasetVar)
    chromosomes <- tclvalue(chromosomesVar)

    if (chromosomes == "") {
      chromosome_list <- NULL
    } else {
      chromosome_list <- strsplit(chromosomes, ", ?")[[1]]
    }

    # Execute the provided code for download
    biomart <- BioMartData$new(biomart = mart, dataset = dataset)
    biomart$get_data(chromosomes = chromosome_list)
   
    # Add some logic here to show the results or log messages if needed
  }

  btn <- tkbutton(tt, text = "Download", command = onDownloadClick)

  # Position controls in the window using grid layout with padding
  tkgrid(lab1, ent1, padx=10, pady=10)
  tkgrid(lab2, ent2, padx=10, pady=10)
  tkgrid(lab3, ent3, padx=10, pady=10)
  tkgrid(btn, padx=10, pady=20)

  tkfocus(tt)
}

############################################################################################

library(tcltk)

biomartDownloadWindow <- function() {
  # Create a new top-level window
  tt <- tktoplevel()
  tkwm.title(tt, "DigestR BioMart Downloader")
  tkwm.geometry(tt, "400x250")  # Set the window size

  # Tcl variables for entry widgets
  biomartVar <- tclVar("genes")
  datasetVar <- tclVar("btaurus_gene_ensembl")
  chromosomesVar <- tclVar("1, 2")
  searchPatternVar <- tclVar("taurus")

  # Create and position the controls
  lab1 <- tklabel(tt, text = "Biomart")
  ent1 <- tkentry(tt, textvariable = biomartVar, width=30)
  
  lab2 <- tklabel(tt, text = "Dataset")
  ent2 <- tkentry(tt, textvariable = datasetVar, width=30)
  
  lab3 <- tklabel(tt, text = "Chromosomes")
  ent3 <- tkentry(tt, textvariable = chromosomesVar, width=30)

  lab4 <- tklabel(tt, text = "Search Pattern")
  ent4 <- tkentry(tt, textvariable = searchPatternVar, width = 30)

  # Function to perform the dataset search
  onSearchClick <- function() {
  search_pattern <- tclvalue(searchPatternVar)
  
  if (search_pattern != "") {
    search_results <- searchDatasets(mart = ensembl, pattern = search_pattern)
    # Display search results to the user (you can decide how to display this information)
    # For example, you can create a new window or use a message box
    tkmessageBox(message = paste("Search Results:\n", paste(search_results, collapse = "\n")))
  }
}

#searchBtn <- tkbutton(tt, text = "Search Datasets", command = onSearchClick)

# Position controls for search pattern
#tkgrid(lab4, ent4, padx = 10, pady = 10)
#tkgrid(searchBtn, padx = 10, pady = 20)

  # When the button is clicked, this function will be executed
  onDownloadClick <- function() {
    mart <- tclvalue(biomartVar)
    dataset <- tclvalue(datasetVar)
    chromosomes <- tclvalue(chromosomesVar)

    if (chromosomes == "") {
      chromosome_list <- NULL
    } else {
      chromosome_list <- strsplit(chromosomes, ", ?")[[1]]
    }

    # Execute the provided code for download
    biomart <- BioMartData$new(biomart = mart, dataset = dataset)
    biomart$get_data(chromosomes = chromosome_list)
   
    # Add some logic here to show the results or log messages if needed
  }

  btn <- tkbutton(tt, text = "Download proteome", command = onDownloadClick)
  searchBtn <- tkbutton(tt, text = "Search Datasets", command = onSearchClick)

  # Position controls in the window using grid layout with padding
  tkgrid(lab1, ent1, padx=10, pady=10)
  tkgrid(lab2, ent2, padx=10, pady=10)
  tkgrid(lab3, ent3, padx=10, pady=10)
  tkgrid(lab4, ent4, padx=10, pady=10)
  tkgrid(btn, padx=10, pady=20)
  tkgrid(searchBtn, padx=20, pady=20)

  tkfocus(tt)
}

#########################################################################
biomartDownloadWindow <- function() {
  # Create a new top-level window
  tt <- tktoplevel()
  tkwm.title(tt, "DigestR BioMart Downloader")
  tkwm.geometry(tt, "400x250")  # Set the window size

  # Tcl variables for entry widgets
  biomartVar <- tclVar("genes")
  datasetVar <- tclVar("btaurus_gene_ensembl")
  chromosomesVar <- tclVar("1, 2")
  searchPatternVar <- tclVar("taurus")

  # Create and position the controls
  lab1 <- tklabel(tt, text = "Biomart")
  ent1 <- tkentry(tt, textvariable = biomartVar, width=30)
  
  lab2 <- tklabel(tt, text = "Dataset")
  ent2 <- tk2combobox(win, values = options)
  
  lab3 <- tklabel(tt, text = "Chromosomes")
  ent3 <- tkentry(tt, textvariable = chromosomesVar, width=30)

  lab4 <- tklabel(tt, text = "Search Datasets")
  ent4 <- tkentry(tt, textvariable = searchPatternVar, width = 30)

  # Function to perform the dataset search
# Function to perform the dataset search
onSearchClick <- function() {
  search_pattern <- tclvalue(searchPatternVar)
  
  if (search_pattern != "") {
    # Initialize the ensembl object and retrieve datasets
    # Add message upon clicking
    ensembl <- useEnsembl(biomart = "genes")
    datasets <- listDatasets(ensembl)
    
    # Perform the dataset search
    search_results <- searchDatasets(mart = ensembl, pattern = search_pattern)
    
    # Display search results to the user
    tkmessageBox(message = paste("Search Results:\n", paste(search_results, collapse = "\n")))
    log_message("Fetching Datasets...")
  }
}

  # Create the search button
  searchBtn <- tkbutton(tt, text = "Search Datasets", command = onSearchClick)

  # When the button is clicked, this function will be executed
  onDownloadClick <- function() {
    mart <- tclvalue(biomartVar)
    dataset <- tclvalue(datasetVar)
    chromosomes <- tclvalue(chromosomesVar)

    if (chromosomes == "") {
      chromosome_list <- NULL
    } else {
      chromosome_list <- strsplit(chromosomes, ", ?")[[1]]
    }

    # Execute the provided code for download
    biomart <- BioMartData$new(biomart = mart, dataset = dataset)
    biomart$get_data(chromosomes = chromosome_list)
   
    # Add some logic here to show the results or log messages if needed
  }

  # Create the download button
  btn <- tkbutton(tt, text = "Download proteome", command = onDownloadClick)

  # Position controls in the window using grid layout with padding
  
tkgrid(lab1, ent1, padx=10, pady=10)
tkgrid(lab2, ent2, padx=10, pady=10)
tkgrid(lab3, ent3, padx=10, pady=10)
tkgrid(lab4, ent4, padx=10, pady=10)

# Create a new frame to hold the buttons
buttonFrame <- tkframe(tt)
tkgrid(buttonFrame, columnspan=2, padx=10, pady=10)

# Position the buttons side by side in the frame
tkgrid(btn, column=1, row=5, padx=10, pady=20, sticky="w")
tkgrid(searchBtn, column=0, row=5, padx=10, pady=20, sticky="e")

  tkfocus(tt)
}

# biomartDownloadWindow()
############################################################################

