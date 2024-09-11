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

#' vd()
#'
#' This function provides a graphical user interface (GUI) for generating Venn diagrams 
#' from peptide sequences stored in CSV files. Users can select multiple CSV files from 
#' a directory, and the Venn diagram will display common and unique peptides across the 
#' selected files. The resulting Venn diagram can be saved as a PDF file.
#'
#' @return A graphical user interface that allows users to upload CSV files, select files 
#' for analysis, and generate a Venn diagram of peptide sequences.
#'
#' @details
#' The function offers the following main features:
#' \itemize{
#'   \item \strong{File Upload:} Users can select a directory containing CSV files and 
#'   choose specific files for Venn diagram generation.
#'   \item \strong{Venn Diagram Creation:} The function generates a Venn diagram based 
#'   on the unique peptide sequences found in each selected CSV file. The diagram displays 
#'   intersections between groups of peptides.
#'   \item \strong{Saving Output:} The resulting Venn diagram is saved as a PDF file. 
#'   The user can choose the file name and location for the saved output.
#' }
#'
#' @section CSV File Format:
#' Each CSV file must contain a column named \code{pep_seq}, which includes the peptide 
#' sequences. Only the unique peptide sequences in each file will be used for generating 
#' the Venn diagram.
#'
#' @section Instructions:
#' \itemize{
#'   \item Upload CSV files from a directory using the "Browse" button.
#'   \item If no files are selected, the function will use all the CSV files in the chosen directory.
#'   \item After selecting the files, click the "Generate Venn Diagram" button to create the plot.
#'   \item The Venn diagram will be saved as a PDF at the location specified by the user.
#' }
#'
#' @section Dependencies:
#' The function depends on the following R packages:
#' \itemize{
#'   \item \code{tcltk}: For creating the graphical user interface.
#'   \item \code{VennDiagram}: For generating Venn diagrams.
#'   \item \code{readr}: For reading CSV files.
#'   \item \code{grid}: For displaying the generated Venn diagram.
#'   \item \code{tkrplot}: For interactive plotting.
#' }
#'
#' @examples
#' \dontrun{
#' venn_diagram_generator()
#' }
#'
#' @export

vd <- venn_diagram_generator

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

#' zm()
#' Interactive GUI for zooming and scrolling
#' Zoom and Manipulate View using zoom
#'
#' This function opens a graphical user interface (GUI) that allows users to manipulate the zoom and view of a plot.
#' Users can choose between zooming or scrolling mode, adjust the zoom increment, and perform various zooming actions,
#' such as zooming in, zooming out, and focusing on specific regions of the plot. The GUI also provides options for
#' resetting the view to full, centering the view, and more.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
zm <- zoom

#' ol()
#' Overlay Spectra and Generate Plot
#'
#' This function overlays spectra based on the provided list of spectrum indices and
#' generates a plot with the overlaid spectra. It supports different types of overlay
#' palettes and allows for customization of the overlay behavior.
#'
#' @param askUsr A logical value indicating whether to prompt the user for overlay list creation.
#'               Default is \code{TRUE}.
#' @param offset A numeric value indicating the offset between overlaid spectra. If not provided,
#'               the offset from global settings is used.
#' @param ... Additional parameters passed to internal plotting functions.
#'
#' @return None (invisible return).
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item It defines default overlay palettes for different numbers of colors.
#'   \item It determines the current spectrum and its dimensions.
#'   \item If \code{askUsr} is \code{TRUE} or \code{overlayList} is \code{NULL}, it opens
#'     the GUI for creating an overlay list and returns.
#'   \item If \code{askUsr} is \code{FALSE} and an overlay list is present, it overlays
#'     the spectra as specified by the overlay list.
#'   \item The function handles plot text overlay, removing the current spectrum, and plotting
#'     the overlaid spectra.
#' }
#'
#' @examples
#' \dontrun{
#'   # Overlay spectra based on the provided list
#'   ol(askUsr = FALSE, offset = 0.1)
#' }
#'
#' @importFrom graphics legend par plot
#' @export
ol <- overlay 

#' ct()
#' Wrapper function, displays the plot settings panes in ps()
#' Plot Settings GUI
#'
#' This function creates a graphical user interface (GUI) window for adjusting plot settings for DIANA spectra.
#' The GUI allows users to change plot colors, configure color options for various components, and switch
#' between different settings panels.
#'
#' @param dispPane A character specifying the default displayed pane in the GUI. Options are "co" (default) for
#'                 the "Plot Colors" pane, and "sp" for the "Spectra" pane.
#'
#' @return None (invisible return).
#'
#' @importFrom tcltk tkwm.title tkfocus tkwm.deiconify tkadd tkgrid ttksizegrip tkselect tkcurselection tkitemconfigure tkselection.set tkbind tkdestroy
#' @export
ct <- plot_settings

#' co()
#' Wrapper function, displays the plot colors pane in ps()
#' @param dispPane A character specifying the default displayed pane in the GUI. Options are "co" (default) for
#'                 the "Plot Colors" pane, and "sp" for the "Spectra" pane.
#'
#' @return None (invisible return).
#'
#' @importFrom tcltk tkwm.title tkfocus tkwm.deiconify tkadd tkgrid ttksizegrip tkselect tkcurselection tkitemconfigure tkselection.set tkbind tkdestroy
#' @export
co <- plot_colors 

#' gl()
#' Gene Name Threshold GUI
#'
#' This function opens a graphical user interface (GUI) that allows users to set threshold values for gene names.
#' Users can interact with the GUI to select files, set and apply threshold values, and reset to default values.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
glab <- gene_labeling

#' dd()
#' Refreshes the main plot without changing settings
#' @import tcltk2
#'
#' @export
dd <- redraw

#' Sa()
#' Save File As DIANA Compressed File
#'
#' This function allows the user to save the current DIANA file in a compressed format.
#'
#' @param saveFileName Character string containing the desired save file name. If not provided,
#'                     a file dialog will prompt the user to choose a save location.
#' @return None (invisible return).
#'
#' @export
sa <- save

#' ud()
#' Undo Last Action
#'
#' This function allows users to undo the most recent action performed in the application. It restores the state of various global variables to their previous values, effectively reversing the most recent changes made.
#'
#' @return NULL (The function primarily works with global variables and their modifications.)
#'
#' @details
#' The \code{undo} function works by decrementing the undo index in the \code{oldFolder} variable, which stores the history of changes. It resets each of the global files and adds new zoom changes to the \code{oldFolder} index.
#' 
#' @importFrom tcltk2
#' @export
undo <- ud
