# Aliases

#' gp()
#' Generate Proteome Data Using DigestR BioMart Downloader
#'
#' This function creates a graphical user interface (GUI) using tcltk2 to allow users
#' to interactively search for and download proteome data from BioMart using the DigestR package.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
gp <- generate_proteome

#'   cs()
#'   Plot the cut sites of a given protease on a protein sequence.
#'
#'   Parameters:
#'   - prot (character): Name of the protease.
#'   - colour (character): Color for the lines on the plot.
#' 
#'   Returns:
#'   - Invisible: Returns the protease name invisibly.
#' 
#'   This function reads protease cut site data from a CSV file and visualizes
#'   the cut sites of the specified protease on a protein sequence. It retrieves
#'   the current gene information and finds the corresponding sequence based on
#'   the gene index. Then, it identifies the cut sites of the specified protease
#'   within the sequence and plots vertical lines at those positions using the
#'   provided color.
#' 
#'   Note: The function assumes that the necessary data structures like
#'   'globalSettings', 'species', and 'protCutSites' are available in the
#'   environment.
#' 
#'   Example usage:
#'   plotCutSite("Trypsin", "black")
#'   """
cs <- display_protease_cut_sites

#' csd()
#'
#' This function generates various plots to visualize cleavage site distributions from
#' input CSV files. The plots include percentages of mean first and last letter counts
#' for different groups of data.
#'
#' @param csv_dir A character string specifying the directory containing the input CSV files.
#' @param plot_type A character string specifying the type of plot to generate. Options are
#'                  "Nter" (N-terminal cleavage sites), "Cter" (C-terminal cleavage sites),
#'                  or "Combined" (combined plots for mean first and last letter counts).
#'
#' @return A ggplot2 object representing the generated plot.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The CSV files in the specified directory are grouped based on a second string
#'     in their filenames.
#'   \item First and last letters are extracted from each sequence in the CSV files.
#'   \item Percentages of mean first and last letter counts are calculated for each group.
#'   \item The function creates plots based on the specified plot type: "Nter" (N-terminal),
#'     "Cter" (C-terminal), or "Combined" (both mean first and last letter counts).
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate a plot for N-terminal cleavage sites
#'   nter_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Nter")
#'
#'   # Generate a plot for C-terminal cleavage sites
#'   cter_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Cter")
#'
#'   # Generate a combined plot for mean first and last letter counts
#'   combined_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Combined")
#' }
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar ggtitle xlab ylab theme_bw
#' @export
csd <- cut_sites_distribution

#' pd()
#' Create Density Plot for Peptide Length Distribution
#'
#' This function generates density plots for the distribution of peptide lengths across
#' different groups based on input CSV files. The function provides options to create
#' different types of plots including overlay density plot, ridges plot, and colored
#' ridges plot.
#'
#' @param csv_dir A character string specifying the directory where the CSV files are located.
#' @param plot_type A character string specifying the type of density plot to generate.
#'                  Options are "overlay", "ridges", or "colored_ridges".
#'
#' @return A ggplot2 object representing the generated density plot.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The CSV files in the specified directory are grouped based on a second string
#'     in their filenames.
#'   \item The function calculates the number of characters in the 'pep_seq' column for
#'     each group and calculates the mean length.
#'   \item Depending on the \code{plot_type} argument, the function creates an overlay
#'     density plot, a ridges plot, or a colored ridges plot using the \pkg{ggplot2}
#'     and \pkg{ggridges} packages.
#'   \item The mean length is annotated on the plot.
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate an overlay density plot
#'   overlay_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "overlay")
#'
#'   # Generate a ridges plot
#'   ridges_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "ridges")
#'
#'   # Generate a colored ridges plot
#'   colored_ridges_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "colored_ridges")
#' }
#'
#' @importFrom ggplot2 ggplot geom_density labs theme_bw geom_vline geom_text scale_fill_viridis_c annotate
#' @importFrom ggridges geom_density_ridges2 geom_density_ridges_gradient
#' @importFrom dplyr group_by summarise
#' @export
pd <- peptides_distribution

#' up()
#' Find and Save Unique Peptides from a given mascot file
#'
#' This function prompts the user to select directories and files, reads input CSV files,
#' and identifies unique peptides present in one file but not in another. It then saves
#' these unique peptides to a new CSV file.
#'
#' @return A character string indicating the path of the destination CSV file
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The user is prompted to select a directory to work in.
#'   \item The user is prompted to select an input CSV file yo use as a Query (L1).
#'   \item The user is prompted to select another input CSV file as the experimental group (L2).
#'   \item The CSV files are read with headers assumed to be present.
#'   \item The function checks if the 'pep_seq' column exists in both L1 and L2 files.
#'     If not, an error is raised.
#'   \item Unique peptides are identified in L2 that are not present in L1.
#'   \item A new filename is created for the output CSV file by appending "_unique" to
#'     the L2 file name (without extension).
#'   \item The unique peptides are saved to the destination CSV file.
#' }
#'
#' @examples
#' \dontrun{
#'   # Call the function
#'   uniquePeptides()
#' }
#'
#' @importFrom utils read.csv write.csv
#' @export
up <- unique_peptides

#' rd()
#' Remove Duplicate Entries from a CSV File
#'
#' This function prompts the user to select an input CSV file, reads its content,
#' removes duplicate entries based on a specified column, and writes the unique
#' data to a destination CSV file.
#'
#' @return A character string indicating the path of the destination CSV file
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The user is prompted to select a directory, which becomes the working directory.
#'   \item The user is prompted to select the input CSV file.
#'   \item The CSV file is read, skipping the first 3 rows, and the header is assumed to be present.
#'   \item The function checks if the 'pep_seq' column exists in the data. If not, an error is raised.
#'   \item Duplicate values in the 'pep_seq' column are removed, and the resulting data frame is stored.
#'   \item The user is prompted to select a destination CSV file.
#'   \item The unique data is written to the destination CSV file, excluding row names.
#' }
#'
#' @examples
#' \dontrun{
#'   # Call the function
#'   RemoveDuplicate()
#' }
#'
#' @importFrom utils read.csv write.csv
#' @export	 
rd <- prepare_mascot

#' pm()
#'
#' This function provides a GUI for processing Mascot files. 
#' The user can select a species and choose to either convert a single file or all files in a folder 
#' and its subfolders.
#'
#' @return Invisible. The function internally updates global settings based on user inputs 
#' from the GUI and then calls a batch conversion function.
#'
#' @note The function uses the `tcltk` package to create the GUI.
#'
#' @examples
#' \dontrun{
#' pm()
#' }
#'
#' @import tcltk
#' @export
pm <- process_mascot

#' fo()
#' Open and Load Spectrum Files for `digestR`
#'
#' The `fo` function is responsible for loading one or more spectrum files, processing them, 
#' and integrating them into the `digestR` environment. If the user does not specify a 
#' filename, a GUI-based file selector will be triggered. The function will then read the 
#' specified files, verify their format and content, and make necessary updates to the global 
#' state. This function also provides error handling, ensuring that files with unsupported or 
#' corrupted formats do not crash the application.
#' 
#' @param fileName A character string specifying the path to a spectrum file to be loaded. 
#' If `fileName` is missing, the user will be prompted to select files via a GUI-based selector.
#' @param ... Additional arguments passed to other methods or functions.
#' 
#' @return Returns an invisibly the user-selected or provided file paths.
#' 
#' @examples 
#' \dontrun{
#' # Load a single spectrum file (assumes valid path and format)
#' fo("path_to_spectrum_file.txt")
#' 
#' # Trigger GUI-based file selector for user to select files
#' fo()
#' }
#' 
#' @seealso `createObj`, `dianaHead`, `log_message`, `myAssign`, `refresh`
#'
#' @export#'
#' The `fo` function is responsible for loading one or more spectrum files, processing them, 
#' and integrating them into the `digestR` environment. If the user does not specify a 
#' filename, a GUI-based file selector will be triggered. The function will then read the 
#' specified files, verify their format and content, and make necessary updates to the global 
#' state. This function also provides error handling, ensuring that files with unsupported or 
#' corrupted formats do not crash the application.
#' 
#' @param fileName A character string specifying the path to a spectrum file to be loaded. 
#' If `fileName` is missing, the user will be prompted to select files via a GUI-based selector.
#' @param ... Additional arguments passed to other methods or functions.
#' 
#' @return Returns an invisibly the user-selected or provided file paths.
#' 
#' @examples 
#' \dontrun{
#' # Load a single spectrum file (assumes valid path and format)
#' fo("path_to_spectrum_file.txt")
#' 
#' # Trigger GUI-based file selector for user to select files
#' fo()
#' }
#' 
#' @seealso `createObj`, `dianaHead`, `log_message`, `myAssign`, `refresh`
#'
#' @export
fo <- file_open

