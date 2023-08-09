file_path <- system.file("logger.R", package = "digestR")
print(file_path)
if (file.exists(file_path)) {
    base::source(file_path)
} else {
    stop("logger.R was not found!")
}

library(shiny)

library(promises)
library(future)
library(ggplot2)

# Set the plan to multicore to allow for parallel execution
plan(multisession)
#plan(sequential)

# Define the directories at the top of your script
proteomesDirectoryPath <- "data/proteomes"
peptidesDirectoryPath <- "data/peptides"
dcfDirectoryPath <- "data/dcf"


maybeCreateDirectory <- function(dir) {
  # Check if the directory exists
  if (!dir.exists(dir)) {
    # If the directory does not exist, create it
    dir.create(dir, recursive = TRUE)  # The recursive = TRUE argument means that it will create any necessary parent directories as well
  }
}


create_dcf_path <- function(peptidesFile) {
  # Get the filename (without the extension)
  filename <- tools::file_path_sans_ext(basename(peptidesFile))
  
  # Add the ".dcf" extension to the filename
  filename_dcf <- paste0(filename, ".dcf")
  
  # Combine the directory path and the new filename to get the full path
  full_path_dcf <- file.path(dcfDirectoryPath, filename_dcf)
  
  return(full_path_dcf)
}



processFile_Soren <- function(proteomeFile, peptidesFile) {
  proteomeFile <- normalizePath(proteomeFile)
  peptidesFile <- normalizePath(peptidesFile)
  saveFileName <- normalizePath(create_dcf_path(peptidesFile), mustWork = FALSE)

  log_message(paste('Preparing proteome file', proteomeFile, sep=' '))


  if (proteomeFile == "" || is.null(proteomeFile)) {
    stop("proteomeFile is missing or empty")
  }
  if (!file.exists(proteomeFile)) {
    stop("proteomeFile does not exist at the provided path")
  }

  if (peptidesFile == "" || is.null(peptidesFile)) {
    stop("peptidesFile is missing or empty")
  }
  if (!file.exists(peptidesFile)) {
    stop("peptidesFile does not exist at the provided path")
  }

  log_message(paste('Preparing peptides file', peptidesFile, sep=' '))
  log_message(paste('Results will be saved in', saveFileName, sep=' '))

  lSpecies <- loadSpecies(proteomeFile)
  log_message("loadSpecies done")

  future({

    print('Checking variables in future() block...')
    print(paste('Peptides file', peptidesFile, sep=' '))
    print(paste('Results will be saved in', saveFileName, sep=' '))

    print('Processing file...')

    # This code will be executed in a separate R process
  
    result <- processFile(peptidesFile, lSpecies)

    log_message(paste('Done processing file: ', peptidesFile))  

    if(is.numeric(result)) {
      if(result == -1) {
        print(paste0('Mascot file could not be read or found, ', saveFileName, " could not be generated."))
      } else {
        print(paste0('No matches found, ', saveFileName, " could not be generated."))
      }
    } else {
      sName <- writeDIANA(saveFileName, lSpecies, result) 
      print(paste0(sName, ' saved at ', format(Sys.time(), "%H:%M"), '.'))
    }
  }) %...>% {
    # This code will be executed in the main R process, after the future has resolved
    # You could put any post-processing code here
  } %...!% {
    # This code will be executed in the main R process if the future throws an error
    error <- .

    log_message(paste("An error occurred while processing the file:", error$message, 
                      "\nTraceback:\n", paste(capture.output(traceback()), collapse = "\n")))
    log_message(paste('Proteome file', proteomeFile, sep=' '))
    log_message(paste('Peptides file', peptidesFile, sep=' '))
  }
}


set_directory <- function(dir) {
  # Check if the directory exists
  if (!dir.exists(dir)) {
    # If the directory does not exist, create it
    dir.create(dir, recursive = TRUE)  # The recursive = TRUE argument means that it will create any necessary parent directories as well
  }

  # Set the working directory
  setwd(dir)

  log_message('Working directory:', getwd())
}

# Define reactive expression for globalSettings
initializeGlobalSettings <- function() {
  reactiveValues(
    proteomeFiles = list(),
    peptidesFiles = list(),
    results = NULL
  )
}


handleFileUpload <- function(input, fileInputName, globalSettings, globalSettingsName, label, uploadDirectory) {
  observe({
    inFile <- input[[fileInputName]]
    
    if (is.null(inFile)) {
      return(NULL)
    }
    
    # Create the upload directory if it doesn't already exist
    if (!dir.exists(uploadDirectory)) {
      dir.create(uploadDirectory, recursive = TRUE)
    }
    
    # Move the uploaded file to the upload directory
    fileDestination <- file.path(uploadDirectory, inFile$name)
    file.copy(from = inFile$datapath, to = fileDestination, overwrite = TRUE)

    # Update globalSettings using reactive conduct
    globalSettings[[globalSettingsName]][[inFile$name]] <- list(file = fileDestination, label = label)
  })
}


createShinyUI <- function() {
  fluidPage(
    titlePanel("digestR"),
    fluidRow(
      column(3,
        wellPanel(
          #style = "padding: 15px; height: 90vh; max-width: 400px; overflow-y: auto; background-color: #f5f5f5;",
          uiOutput("proteome_selector"),
          uiOutput("peptides_selector"),
          actionButton("button_run", "Run"),

        )
      ),
      column(9,
        tabsetPanel(

          tabPanel("Results",
            # Contents for displaying results...
            h4("Log Messages:"),
            verbatimTextOutput("log"),                
          ),

          tabPanel("Upload & BioMart",
            fluidRow(
              column(12,
                wellPanel(
                  tags$h2("Upload"),  # Add a title for the panel

                  fileInput("file_proteome", "Choose Proteome CSV File to upload",
                            multiple = TRUE,
                            accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                  fileInput("file_peptides", "Choose Peptides CSV File to upload",
                            multiple = TRUE,
                            accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                ),
                wellPanel(
                  tags$h2("BioMart"),  # Add a title for the panel
                  textInput("text_biomart", label = "Biomart", value = "ensembl"),
                  textInput("text_biomart_dataset", label = "Dataset", value = "btaurus_gene_ensembl"),
                  textInput("text_biomart_chromosomes", label = "Chromosomes", value = "1, 2"),
                  actionButton("button_download", "Download"),
                  verbatimTextOutput("settingsOutput")
                )
              )
            )
          ),



          tabPanel("Plotting",
            # Contents for displaying results...
            selectInput("var", 
                  label = "Choose a variable", 
                  choices = colnames(mtcars), 
                  selected = "mpg"),
            plotOutput("densityPlot")
          ),


          tabPanel("Global Settings",
            # Contents for displaying global settings...
            verbatimTextOutput("globalSettingsOutput"),
        
          ),
          
        )
      )
    )
  )
}


start_app <- function() {
  # Initialize globalSettings
  set_directory("~/digestR")

  maybeCreateDirectory(dcfDirectoryPath)
  maybeCreateDirectory(proteomesDirectoryPath)
  maybeCreateDirectory(peptidesDirectoryPath)

  globalSettings <- initializeGlobalSettings()
  
  ui <- createShinyUI()

  server <- function(input, output, session) {
    options(shiny.maxRequestSize = 500*1024^2)  # Set upload limit to 500MB
      # Handle file uploads
    handleFileUpload(input, "file_proteome", globalSettings, "proteomeFiles", "Choose Proteome CSV File", proteomesDirectoryPath)
    handleFileUpload(input, "file_peptides", globalSettings, "peptidesFiles", "Choose Peptides CSV File", peptidesDirectoryPath)

    output$workingDirectory <- renderText({
      getwd()
    })   

    output$settingsOutput <- renderPrint({
      # Initialize an empty list to store settings
      settings <- list()

      # If any species files have been uploaded, add their labels to the settings list
      if (length(globalSettings$proteomeFiles) > 0) {
        for (i in seq_along(globalSettings$proteomeFiles)) {
          settings[paste0("Species File ", i)] <- names(globalSettings$proteomeFiles)[[i]]
        }
      }

      # If any peptide files have been uploaded, add their labels to the settings list
      if (length(globalSettings$peptidesFiles) > 0) {
        for (i in seq_along(globalSettings$peptidesFiles)) {
          settings[paste0("Peptide File ", i)] <- names(globalSettings$peptidesFiles)[[i]]
        }
      }

      # If no files have been uploaded, display a message indicating this
      if (length(settings) == 0) {
        settings["No files"] <- "No files have been uploaded yet."
      }

      settings['results'] <- class(globalSettings$results)

      # Print the formatted output
      cat(paste(names(settings), settings, sep = " = ", collapse = "\n"))
    })

    observeEvent(input$button_download, {
      log_message('Downloading proteome')

      mart = input$text_biomart
      dataset = input$text_biomart_dataset
      chromosomes = input$text_biomart_chromosomes

      if (chromosomes == "") {
        chromosome_list <- NULL
      } else {
        chromosome_list <- strsplit(chromosomes, ", ?")[[1]]
      }

      biomart <- BioMartData$new(biomart = "ensembl", dataset = "btaurus_gene_ensembl")

      # Retrieve and process the data
      biomart$get_data(chromosomes = chromosome_list)
    }) 

    # Define reactive expressions for the directories
    proteomesDirectory <- reactivePoll(1000, session, 
      checkFunc = function() {
        # This function is called to check if the value has changed
        # We use file.info to get information about the directory, 
        # including the last modification time
        file.info(proteomesDirectoryPath)$mtime
      }, 
      valueFunc = function() {
        # This function is called to get the actual value (the list of files)
        list.files(proteomesDirectoryPath)
      }
    )

    peptidesDirectory <- reactivePoll(1000, session, 
      checkFunc = function() {
        file.info(peptidesDirectoryPath)$mtime
      }, 
      valueFunc = function() {
        list.files(peptidesDirectoryPath)
      }
    )

    output$proteome_selector <- renderUI({
      selectInput("proteome_selector", "Select proteome files", choices = proteomesDirectory(), multiple = TRUE)
    })

    output$peptides_selector <- renderUI({
      selectInput("peptides_selector", "Select peptides files", choices = peptidesDirectory(), multiple = TRUE)
    })

    # Then, use these variables in your code
    observeEvent(input$proteome_selector, {
      # Update the proteomeFiles element in globalSettings
      selected_files <- input$proteome_selector
      full_paths <- file.path(proteomesDirectoryPath, selected_files)
      globalSettings$proteomeFiles <- full_paths
    })

    observeEvent(input$peptides_selector, {
      # Update the peptideFiles element in globalSettings
      selected_files <- input$peptides_selector
      full_paths <- file.path(peptidesDirectoryPath, selected_files)
      globalSettings$peptidesFiles <- full_paths
    })

    observeEvent(input$button_run, {  # This code gets executed when the button is clicked

      proteomeFile = globalSettings$proteomeFiles[[1]]
      peptidesFiles = globalSettings$peptidesFiles
      
      msg <- paste("Selected proteomeFile:", proteomeFile)
      log_message(msg)
      
      msg <- paste("Selected peptidesFile:", peptidesFiles)
      log_message(msg)
     
      #showModal(modalDialog(
      #  title = "Please wait",
      #  "Processing...",
      #))


      for (peptidesFile in peptidesFiles) {
        log_message(paste('Processing', peptidesFile, sep=' '))
        processFile_Soren(proteomeFile, peptidesFile)
      }

      #removeModal()

    })

    output$densityPlot <- renderPlot({
      ggplot(mtcars, aes_string(x = input$var)) +
        geom_density(fill = 'blue', alpha = 0.5) +
        theme_minimal() +
        labs(x = input$var, y = "Density", title = paste("Density Plot of", input$var))
    })

    output$globalSettingsOutput <- renderPrint({
      # Initialize an empty list to store settings
      settings <- list()
      
      # For each element in globalSettings, add it to the settings list
      for (name in names(globalSettings)) {
        value <- globalSettings[[name]]
        settings[[name]] <- value
      }
      
      # Print the formatted output
      cat(paste(names(settings), sapply(settings, function(x) paste(capture.output(print(x)), collapse = "\n")), sep = " = ", collapse = "\n"))
    })

    # Read the log file and send it to the UI
    output$log <- renderText({
      # reactivePoll checks every 1000 milliseconds (1 second)
      logFileContent <- reactivePoll(1000, session,
        checkFunc = function() {
          # The check function returns the last modification time of the log file
          if (file.exists("digestR.log")) {
            x <- file.info("digestR.log")$mtime
            return(x)
          }
        },
        valueFunc = function() {
          # The value function returns the contents of the log file
          if (file.exists("digestR.log")) {
            lines <- readLines("digestR.log")
            # Reverse the order of lines and only take the last 15
            lines <- tail(lines, 80)
            paste(lines, collapse = "\n")
          } else {
            "No log messages."
          }
        }
      )

      logFileContent()
    })

  }
  shinyApp(ui = ui, server = server)
}

# Run the Shiny app
#start_app()
