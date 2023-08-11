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

