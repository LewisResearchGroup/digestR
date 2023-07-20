# Define reactive expression for globalSettings
initializeGlobalSettings <- function() {
  reactiveValues(
    speciesFiles = list(),
    peptideFiles = list(),
    results = NULL
  )
}

handleFileUpload <- function(input, fileInputName, globalSettings, globalSettingsName, label) {
  observe({
    inFile <- input[[fileInputName]]
    
    if (is.null(inFile)) {
      return(NULL)
    }
    
    # Update globalSettings using reactive conduct
    globalSettings[[globalSettingsName]][[inFile$name]] <- list(file = inFile$datapath, label = label)

    print(globalSettings)

  })
}


createShinyUI <- function() {
  fluidPage(
    titlePanel("digestR"),
    fluidRow(
      column(4, id = "sidebar", style = "padding: 15px; height: 90vh; max-width: 400px; overflow-y: auto; background-color: #f5f5f5;",
        fileInput("file_genome", "Choose Genome CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        fileInput("file_peptides", "Choose Peptides CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        actionButton("button_run", "Run"),  # This creates a button
        verbatimTextOutput("settingsOutput")
      ),
      column(8, 
        textOutput("runOutput")  # This is where the message after clicking Run button will be displayed
      )
    )
  )
}




start_app <- function() {
  # Initialize globalSettings
  globalSettings <- initializeGlobalSettings()
  
  ui <- createShinyUI()

  server <- function(input, output) {
    options(shiny.maxRequestSize = 500*1024^2)  # Set upload limit to 500MB

    # Handle file uploads
    handleFileUpload(input, "file_genome", globalSettings, "speciesFiles", "Choose Genome CSV File")
    handleFileUpload(input, "file_peptides", globalSettings, "peptideFiles", "Choose Peptides CSV File")

    output$settingsOutput <- renderPrint({
      # Initialize an empty list to store settings
      settings <- list()

      # If any species files have been uploaded, add their labels to the settings list
      if (length(globalSettings$speciesFiles) > 0) {
        for (i in seq_along(globalSettings$speciesFiles)) {
          settings[paste0("Species File ", i)] <- names(globalSettings$speciesFiles)[[i]]
        }
      }

      # If any peptide files have been uploaded, add their labels to the settings list
      if (length(globalSettings$peptideFiles) > 0) {
        for (i in seq_along(globalSettings$peptideFiles)) {
          settings[paste0("Peptide File ", i)] <- names(globalSettings$peptideFiles)[[i]]
        }
      }

      # If no files have been uploaded, display a message indicating this
      if (length(settings) == 0) {
        settings["No files"] <- "No files have been uploaded yet."
      }

      # Print the formatted output
      cat(paste(names(settings), settings, sep = " = ", collapse = "\n"))
    })

    observeEvent(input$button_run, {  # This code gets executed when the button is clicked
      output$runOutput <- renderText({
        "Analysing peptides!"
      })

      if (length(globalSettings$speciesFiles) > 0) {
        speciesFile = globalSettings$speciesFiles[[1]]$file
      }

      if (length(globalSettings$peptideFiles) > 0) {
        peptideFile = globalSettings$peptideFiles[[1]]$file
      }

      print('Loading genome.')
      lSpecies <- loadSpecies(speciesFile)
      
      print('Processing file.')
      globalSettings$results <- processFile(peptideFile, lSpecies)

      print('Done')
      print(globalSettings$results)
      print(class(globalSettings$results))
    })
  }

  shinyApp(ui = ui, server = server)
}

# Run the Shiny app
#start_app()
