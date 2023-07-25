# DigestR

**Author**: Dimitri Desmonts de Lamache

## Overview
DigestR is an open-source software developed for the R statistical language environment. The DigestR package is based on rNMR and contains a collection of tools for visualizing and analyzing protein digestion. 
Users can interact interact with DigestR in two major ways: via point and click graphical user interfaces (GUIs) or by entering command directly in the R console. 
This guide is intended to give an overview of DigestR functions.


## Prerequisites

	install.package('biomaRt')
	library(biomarRt)



## Run the following lines in the R console

	library(Rcpp)
	library(tcltk)
	setwd('C:/Users/dimit/OneDrive/Bureau/DigestR')
 
	source('2017_09_16_rNMR_Travis_DD_Edits_Final.R')
	source('digestR_Code_Dimitri_final.R')
	
	

	rm(list = ls())
	devtools::load_all("C:\\Users\\soere\\workspace\\digestR")


## Supported file formats 

DigestR supports Mascot (.csv) and MaxQuant (.mzXML) generated files. These files may be converted to .dcf files using DigestR file conversion functions pm() and im() respectively. 

## GUI Functions Documentation

This document provides details of several Graphical User Interface (GUI) functions implemented in the R programming language using the Tcl/Tk toolkit.

### Converting files. 


### 1. Remove Duplicates Function
The `rd()` function creates a main GUI window, which provides an option to remove duplicate entries from a CSV file. This function first creates a "Remove Duplicates" button. When clicked, the button opens a file dialog for the user to choose a CSV file. After the file is selected, the removeDuplicates function is called to remove duplicate entries based on the 'pep_seq' column. The resulting file is then saved to a location specified by the user.

### 2. Unique Peptides Function
The `up()` function creates a main GUI window, which provides options to select two CSV files (Query and Experimental) and find unique peptides in the Experimental file that are not present in the Query file. The function then writes the unique peptides to a CSV file named with the Experimental filename appended with "_Unique".

### 3. Density Plot Function
The `pd()` function creates a GUI that allows users to create density plots from grouped CSV files located in a specified directory. Users can choose from three types of plots: Overlay, Ridges, and Colored Ridges.

### 4. Cut Sites Distribution
The `csp()` function creates a GUI that allows users to create bar plots from grouped CSV files located in a specified directory. Users can choose from three types of plots: Nter, Cter, and Combined.

### 5. Plotting Cut Sites Function
The `cs()` function creates a GUI window that allows users to select a protease from a dropdown menu and input a color for the lines in the plot. Upon clicking the "Plot" button, the function calls the plotCutSite function, passing the selected protease and color as arguments, to generate a plot of cut sites.

### 6. Testing the Functions
For each function, to test the GUI, call the function without any arguments. For example, to test the up function, simply call up().









---
---
---
---


## Notes



```
### Call the GUI handler function to start the application

    cs()


# rd <- function() {
#   ## creates main window
#   tclCheck()
#   dlg <- tktoplevel()
#   tkwm.title(dlg, 'Prepare Mascot file for digestR')
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
#   
#   ## create Remove Duplicates button
#   onRemoveDuplicates <- function() {
#     source_file <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     # Call the RemoveDuplicate() function passing the selected file
#     result_file <- removeDuplicates(source_file)
#     # Perform actions with the result_file, such as displaying or opening it
#   }
#   removeDuplicatesButton <- ttkbutton(dlg, text='Remove Duplicates', width=15, command=onRemoveDuplicates)
#   
#   ## add button to the main window
#   tkgrid(removeDuplicatesButton, column=1, row=1, padx=5, pady=5)
#   
#   ## make main window stretch when resized
#   tkgrid.columnconfigure(dlg, 1, weight=1)
#   tkgrid.rowconfigure(dlg, 1, weight=1)
#   
#   tkwait.window(dlg)
# }
# 
# removeDuplicates <- function(input_file) {
#   # Read the input CSV file
#   data <- read.csv(input_file, skip = 3, header = TRUE)
#   
#   # Check if 'pep_seq' column exists
#   if (!("pep_seq" %in% colnames(data))) {
#     stop("Error: 'pep_seq' column not found in the input file.")
#   }
#   
#   # Remove duplicate values in the pep_seq column
#   rd_data <- data[!duplicated(data$pep_seq), ]
#   
#   # Prompt the user to select the destination CSV file
#   dest_file <- tclvalue(tcl("tk_getSaveFile", "-initialdir", ".", "-defaultextension", ".csv", "-filetypes", "{{CSV Files} {.csv}}"))
#   
#   # Write the unique data to the destination CSV file
#   write.csv(rd_data, file = dest_file, row.names = FALSE)
#   
#   # Return the path of the destination CSV file
#   return(dest_file)
# }
# 
# # Test the RemoveDuplicate() function
# #rd()
# 
# ################################################################################
# # Unique peptides function  GUI 
# ## Change code so that it uses the current directory to export the file to.
# 
# up <- function() {
#   ## creates main window
#   tclCheck()
#   dlg <- tktoplevel()
#   tkwm.title(dlg, 'Unique peptides')
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
#   
#   ## create select Query button
#   L1Var <- tclVar()  # Variable to store the path of the selected L1 file
#   onSelectL1 <- function() {
#     selectedL1 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     tclvalue(L1Var) <- selectedL1
#   }
#   selectL1Button <- ttkbutton(dlg, text='Query', width=15, command=onSelectL1)
#   
#   ## create select L2 button
#   L2Var <- tclVar()  # Variable to store the path of the selected L2 file
#   onSelectL2 <- function() {
#     selectedL2 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     tclvalue(L2Var) <- selectedL2
#   }
#   selectL2Button <- ttkbutton(dlg, text='Experimental', width=15, command=onSelectL2)
#   
#   ## create Remove Duplicates button
#   onuniquePeptides <- function() {
#     L1 <- read.csv(tclvalue(L1Var), header=TRUE)
#     L2_file <- tclvalue(L2Var)
#     L2 <- read.csv(L2_file, header=TRUE)
#     
#     # Check if 'pep_seq' column exists in L1 and L2
#     if (!("pep_seq" %in% colnames(L1)) || !("pep_seq" %in% colnames(L2))) {
#       stop("Error: 'pep_seq' column not found in the input files.")
#     }
#     
#     # Find unique peptides in L2 not present in L1
#     L2_uniques <- L2[!(L2$pep_seq %in% L1$pep_seq), ]
#     
#     # Get the name of the L2 file without extension
#     l2_filename <- tools::file_path_sans_ext(basename(L2_file))
#     
#     # Create the new filename with "_unique" attached
#     dest_file <- paste0(l2_filename, "_Unique.csv")
#     
#     # Write the unique data to the destination CSV file
#     write.csv(L2_uniques, file = dest_file, row.names = FALSE)
#     
#     tkmessageBox(message=paste("Unique peptides saved to:", dest_file))
#   }
#   uniquePeptidesButton <- ttkbutton(dlg, text='Experimental unique', width=20, command=onuniquePeptides)
#   
#   ## add widgets to the main window
#   tkgrid(selectL1Button, column=1, row=1, padx=5, pady=5)
#   tkgrid(selectL2Button, column=1, row=2, padx=5, pady=5)
#   tkgrid(uniquePeptidesButton, column=1, row=3, padx=5, pady=5)
#   
#   ## make the main window stretch when resized
#   tkgrid.columnconfigure(dlg, 1, weight=1)
#   tkgrid.rowconfigure(dlg, 1, weight=1)
#   
#   tkwait.window(dlg)
# }
# 
# # Test the up() function
# #up()
# 
# ################################################################################
# # Density plot function GUI 
# pd <- function() {
#   # Define variables
#   csv_dir <- ""
#   plot_type <- ""
#   
#   # Function to select the CSV directory
#   selectCSVDir <- function() {
#     csv_dir <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_dir <- tclVar(csv_dir)  # Update the csv_dir variable
#     tkconfigure(csvDirEntry, text = csv_dir)  # Update the csvDirEntry widget
#   }
#   
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     plot_type <<- value
#   }
#   
#   # Function to create the density plot
#   createDensityPlot <- function() {
#     grouped_df <- groupCSVFiles(csv_dir)
#     
#     # Create an empty data frame to store the nchar values for all groups
#     nchar_df <- data.frame(group_name = character(), nchar = integer())
#     
#     for (i in 1:length(grouped_df)) {
#       group_name <- names(grouped_df)[i]
#       pep_seq_col <- grouped_df[[i]]$pep_seq
#       group_nchar_df <- data.frame(group_name = group_name, nchar = nchar(pep_seq_col))
#       nchar_df <- rbind(nchar_df, group_nchar_df)
#     }
#     
#     nchar_mean_df <- nchar_df %>%
#       group_by(group_name) %>%
#       summarise(mean_nchar = mean(nchar))
#     
#     if (plot_type == "Overlay") {
#       plot <- ggplot(nchar_df, aes(x = nchar, group = group_name, fill = group_name)) +
#         geom_density(alpha = 0.5) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         theme_bw() +
#         geom_vline(data = nchar_mean_df, aes(xintercept = mean_nchar, color = group_name), size = 1)
#     } else if (plot_type == "Ridges") {
#       plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = group_name)) +
#         geom_density_ridges2(alpha = 0.5, rel_min_height = 0.01, scale = 7, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         theme_bw()
#     } else if (plot_type == "Colored Ridges") {
#       plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = stat(x))) +
#         geom_density_ridges_gradient(alpha = 0.5, rel_min_height = 0.01, scale = 3, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         scale_fill_viridis_c(name = "Length", option = "H")
#     } else {
#       stop("Invalid plot_type argument. Please choose either 'Overlay', 'Ridges', or 'Colored Ridges'.")
#     }
#     
#     print(plot)
#   }
#   
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Density Plot Creator")
#   
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
#   
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_dir)
#   tkgrid(csvDirEntry, column = 2, row = 1)
#   
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
#   
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
#   
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
#   
#   # Create and configure the Overlay radiobutton
#   overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Overlay", variable = plot_type, value = "Overlay", command = function() setPlotType("Overlay"))
#   tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the Ridges radiobutton
#   ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Ridges", variable = plot_type, value = "Ridges", command = function() setPlotType("Ridges"))
#   tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
#   
#   # Create and configure the Colored Ridges radiobutton
#   coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Colored Ridges", variable = plot_type, value = "Colored Ridges", command = function() setPlotType("Colored Ridges"))
#   tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
#   
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = createDensityPlot)
#   tkgrid(createPlotButton, pady = 10)
#   
#   # Start the event loop
#   tkwait.visibility(win)
#   
# }
# 
# #call the function
# #pd()
# 
# ################################################################################ 
# # Cter, Nter function GUI 
# csp <- function () {
#   
#   # Define variables
#   csv_directory <- ""
#   PlotType <- ""
#   
#   # Function to select the CSV directory
#   selectCSVDir <- function() {
#     csv_directory <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_directory <- tclVar(csv_directory)  # Update the csv_directory variable
#     tkconfigure(csvDirEntry, text = csv_directory)  # Update the csvDirEntry widget
#   }
#   
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     PlotType <<- value
#   }
#   
#   # Function to create the density plot
#   createBarPlot <- function() 
#     #plotCutSites <- function(csv_directory, PlotType) 
#   {
#     library(dplyr)
#     library(ggplot2)
#     
#     extract_first_last_letters_from_csv <- function(file_path) {
#       # Read the CSV file skipping the first 3 lines (assuming headers are in line 4)
#       data <- read.csv(file_path, skip = 0, header = TRUE)
#       
#       # Check if the 'pep_seq' column exists
#       if (!("pep_seq" %in% colnames(data))) {
#         stop("Error: 'pep_seq' column not found in the CSV file.")
#       }
#       
#       # Extract the first and last letters from each sequence
#       data$first_letter <- substr(data$pep_seq, 1, 1)
#       data$last_letter <- substr(data$pep_seq, nchar(data$pep_seq), nchar(data$pep_seq))
#       
#       # Count the occurrences of each first letter
#       first_letter_counts <- table(factor(data$first_letter, levels = c("A", "C", "D", "E", "F", "G","H" ,"K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
#       
#       # Count the occurrences of each last letter
#       last_letter_counts <- table(factor(data$last_letter, levels = c("A", "C", "D", "E", "F", "G","H", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
#       
#       # Combine the counts for first and last letters
#       combined_counts <- merge(
#         data.frame(letter = names(first_letter_counts), first_letter_count = as.numeric(first_letter_counts)),
#         data.frame(letter = names(last_letter_counts), last_letter_count = as.numeric(last_letter_counts)),
#         by = "letter", all = TRUE
#       )
#       combined_counts$combined_count <- rowSums(combined_counts[, c("first_letter_count", "last_letter_count")], na.rm = TRUE)
#       
#       # Replace NA values with 0
#       combined_counts[is.na(combined_counts)] <- 0
#       
#       # Return the combined counts
#       return(combined_counts)
#     }
#     
#     group_csv_files <- function(csv_directory) {
#       # Get the list of CSV files in the directory
#       csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
#       
#       # Create an empty list to store the grouped CSV files
#       grouped_csv <- list()
#       
#       # Loop through the CSV files and group them based on the second string
#       for (i in 1:length(csv_files)) {
#         # Get the filename without the directory path
#         filename <- basename(csv_files[i])
#         # Get the second string in the filename by splitting on "_"
#         second_str <- strsplit(filename, "_")[[1]][2]
#         # If the group already exists in the list, append the CSV file
#         if (second_str %in% names(grouped_csv)) {
#           grouped_csv[[second_str]] <- c(grouped_csv[[second_str]], csv_files[i])
#         }
#         # If the group doesn't exist in the list, create a new list element
#         else {
#           grouped_csv[[second_str]] <- list(csv_files[i])
#         }
#       }
#       
#       # Return the list of grouped CSV files
#       return(grouped_csv)
#     }
#     
#     # Group the CSV files
#     csv <- group_csv_files(csv_directory)
#     
#     # Initialize an empty dataframe to store the results
#     dat <- data.frame()
#     
#     # Iterate over each group of CSV files
#     for (groups in csv) {
#       
#       # Create an empty dataframe to store the results for the group
#       group_results_df <- data.frame(
#         Group = character(),
#         Mean_First_Letter_Count = double(),
#         SD_First_Letter_Count = double(),
#         Mean_Last_Letter_Count = double(),
#         SD_Last_Letter_Count = double(),
#         Mean_Combined_Count = double(),
#         Percentage = double(),
#         SD_Combined_Count = double(),
#         stringsAsFactors = FALSE
#       )
#       
#       # Iterate over each CSV file in the group
#       for (i in groups) {
#         print(i)
#         
#         # Extract the first and last letters from the CSV file
#         file_results <- extract_first_last_letters_from_csv(i)
#         
#         # Get the filename without the directory path
#         filename <- basename(i)
#         
#         # Get the second string in the filename by splitting on "_"
#         grp_str <- rep(strsplit(filename, "_")[[1]][2], nrow(file_results))
#         
#         # Store the label for the group
#         file_results$group <- data.frame(grp_str)
#         
#         # Append the results to the dataframe
#         dat <- rbind(dat, file_results)
#       }
#     }
#     
#     # Aggregate the results by letter and group
#     agg_df <- dat %>%
#       group_by(letter, group) %>%
#       dplyr::summarize(
#         mean_combined_count = mean(combined_count),
#         sd_combined_count = sd(combined_count),
#         se_combined_count = sd_combined_count / sqrt(n()),
#         sum_combined_count = sum(combined_count),
#         mean_first_letter_count = mean(first_letter_count),
#         sd_first_letter_count = sd(first_letter_count),
#         se_first_letter_count = sd_first_letter_count / sqrt(n()),
#         sum_first_letter_count = sum(first_letter_count),
#         mean_last_letter_count = mean(last_letter_count),
#         sd_last_letter_count = sd(last_letter_count),
#         se_last_letter_count = sd_last_letter_count / sqrt(n()),
#         sum_last_letter_count = sum(last_letter_count)
#       )
#     
#     # Calculate percentage of mean_Last_letter_count and mean_first_letter_count with SE
#     agg_df_sum <- agg_df %>%
#       group_by(group) %>%
#       dplyr::summarize(
#         total_sum_combined_count = sum(sum_combined_count),
#         total_sum_first_letter_count = sum(sum_first_letter_count),
#         total_sum_last_letter_count = sum(sum_last_letter_count)
#       )
#     
#     agg_df <- agg_df %>%
#       left_join(agg_df_sum, by = "group") %>%
#       dplyr::mutate(
#         perc_mean_first_letter_count = mean_first_letter_count / total_sum_first_letter_count * 100,
#         perc_mean_last_letter_count = mean_last_letter_count / total_sum_last_letter_count * 100,
#         SE_perc_mean_first_letter_count = sqrt(perc_mean_first_letter_count / 100 * (100 - perc_mean_first_letter_count) / total_sum_first_letter_count),
#         SE_perc_mean_last_letter_count = sqrt(perc_mean_last_letter_count / 100 * (100 - perc_mean_last_letter_count) / total_sum_last_letter_count)
#       )
#     
#     # Convert the grouped dataframe to a regular dataframe
#     agg_df <- data.frame(agg_df)
#     
#     # Convert the group variable to a factor
#     grp <- agg_df$group[[1]]
#     grp <- as.factor(grp)
#     agg_df$group <- grp
#     
#     # Create the appropriate density plot based on the PlotType argument
#     
#     # Nter Percentage Plot
#     if (PlotType == "Nter") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
#               ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "N-terminal cleavage sites",
#              x = "", y = "Percentage") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
#     
#     # Cter Percentage Plot
#     else if (PlotType == "Cter") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
#               ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "C-terminal cleavage sites",
#              x = "", y = "Percentage") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
#     
#     # Plot mean_first_letter_count and mean_last_letter_count as percentage with SE
#     else if (PlotType == "Combined") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
#               ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
#           width = 0.2, position = position_dodge(0.9)
#         ) +
#         ggtitle("Percentage of Mean First Letter Count by Group") +
#         xlab("Letter") +
#         ylab("Percentage of Mean First Letter Count") +
#         theme_bw()
#       
#       ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
#               ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
#           width = 0.2, position = position_dodge(0.9)
#         ) +
#         ggtitle("Percentage of Mean Last Letter Count by Group") +
#         xlab("Letter") +
#         ylab("Percentage of Mean Last Letter Count") +
#         theme_bw()
#     }
#     
#     # If the PlotType argument is not valid
#     else {
#       stop("Error: Invalid PlotType argument. Valid options are 'Nter', 'Cter', and 'Combined'.")
#     }
#     
#     # Return the plot
#     return(cutsiteplots)
#   }
#   
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Cut sites distribution")
#   
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
#   
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_directory)
#   tkgrid(csvDirEntry, column = 2, row = 1)
#   
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
#   
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
#   
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
#   
#   # Create and configure the Overlay radiobutton
#   overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Cter", variable = PlotType, value = "Cter", command = function() setPlotType("Cter"))
#   tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the Ridges radiobutton
#   ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Nter", variable = PlotType, value = "Nter", command = function() setPlotType("Nter"))
#   tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
#   
#   # Create and configure the Colored Ridges radiobutton
#   coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Combined", variable = PlotType, value = "Combined", command = function() setPlotType("Combined"))
#   tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
#   
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = createBarPlot)
#   tkgrid(createPlotButton, pady = 10)
#   
#   # Start the event loop
#   tkwait.visibility(win)
#   
# }
# 
# # call the function
# #csp()
# 
# ################################################################################
# # Plotting cut sites GUI.
# 
# cs <- function() {
#   # Read protease data from CSV
#   protCutSites <- read.csv('C:/Users/dimit/OneDrive/Bureau/digestR/Proteasecutsiteslist.csv')
#   
#   # Create the main window
#   mainWindow <- tktoplevel()
#   tkwm.title(mainWindow, "Plot Cut Site")
#   
#   # Create protease label and dropdown menu
#   proteaseLabel <- tklabel(mainWindow, text = "Select Protease:")
#   tkgrid(proteaseLabel, padx = 10, pady = 10)
#   
#   proteaseDropdown <- tklistbox(mainWindow)
#   tkgrid(proteaseDropdown, padx = 10, pady = 10)
#   
#   # Populate protease dropdown menu with protease names
#   proteaseList <- unique(protCutSites$protease)
#   for (prot in proteaseList) {
#     tkinsert(proteaseDropdown, "end", prot)
#   }
#   
#   # Create color label and entry field
#   colorLabel <- tklabel(mainWindow, text = "Enter Color for Lines:")
#   tkgrid(colorLabel, padx = 10, pady = 10)
#   
#   colorEntry <- tkentry(mainWindow)
#   tkgrid(colorEntry, padx = 10, pady = 10)
#   
#   # Create plot button
#   plotButton <- tkbutton(mainWindow, text = "Plot", command = function() {
#     # Get the selected protease and color
#     selectedProteaseIndex <- tclvalue(tkcurselection(proteaseDropdown))
#     selectedProtease <- tkget(proteaseDropdown, selectedProteaseIndex)
#     selectedColor <- tkget(colorEntry)
#     
#     # Call the plotCutSite function with the selected protease and color
#     plotCutSite(selectedProtease, selectedColor, protCutSites, colorEntry)
#     
#     # Close the GUI window
#     tkdestroy(mainWindow)
#   })
#   tkgrid(plotButton, padx = 10, pady = 10)
# }
# 
# # Function to plot cut sites
# plotCutSite <- function(prot, colour, protCutSites, colorEntry) {
#   # Retrieve current gene information
#   currGene <- globalSettings$geneDisp
#   
#   # Find the index of the current gene in the species data
#   idx <- which(species$genes$name == currGene)
#   
#   if (length(idx) > 0) {
#     # Retrieve the corresponding sequence based on the gene index
#     start <- species$genes$seqStartIdx[idx]
#     end <- start + species$genes$seqLength[idx]
#     currSeq <- substr(species$seq, start, end)
#   } else {
#     stop("No current gene information found. Function doesn't work with proteome-wide view. Call a specific protein.")
#   }
#   
#   # Set gene sequence and colour variables
#   gene <- globalSettings$geneDisp
#   seq <- currSeq
#   colour <- tclvalue(tkget(colorEntry))
#   
#   # Find the index of the protease in the protCutSites data
#   protIndex <- which(as.character(protCutSites$protease) == as.character(prot))
#   col_num <- ncol(protCutSites)
#   
#   # Get the cut sites for the protease
#   protase_cutsite <- protCutSites[protIndex, 3:col_num]
#   loc <- 0
#   
#   # Iterate through the cut sites
#   for (j in protase_cutsite) {
#     if (j == "") {
#       break
#     }
#     print('################################################################')
#     print(j)
#     print('################################################################')
#     
#     prev_index <- 0
#     protcuts <- strsplit(seq, j)
#     
#     # Iterate through the sequence segments
#     for (k in protcuts[[1]]) {
#       k <- paste(k, j, sep = '')
#       print(k)
#       
#       # Find the location of the segment in the sequence
#       loc <- gregexpr(k, seq, fixed = TRUE)
#       start <- loc[[1]][1]
#       len <- nchar(k)
#       
#       # Check if the position is greater than 0 and the segment length is greater than 0
#       if (start > 0 && len > 0) {
#         # Calculate the plot site
#         plot_site <- start + len - 1
#         
#         if (k == "K") {
#           # Check if it is the last character in the sequence
#           if (plot_site == nchar(seq)) {
#             plot_site <- prev_index
#           } else {
#             plot_site <- prev_index + 1
#           }
#           print('IS IT HERE?')
#         }
#         
#         print(plot_site)
#         
#         # Assuming you have initialized a plot or set up the plotting device
#         abline(v = plot_site, col = colour, lty = 2)
#         
#         
#         prev_index <- plot_site
#       }
#     }
#   }
#   return(invisible(prot))
# }
# 
# # Call the GUI handler function to start the application
# #cs()

```
