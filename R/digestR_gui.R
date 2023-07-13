library(tcltk)

#' A Function from gui.R
#'
#' This function does something even more interesting.
#'
#' @export
digestR_gui <- function() {
  
  # Create a new window
  win <- tktoplevel()

  # Create entry fields for the species and file names
  speciesList_entry <- tkentry(win, width = 50)
  speciesFiles_entry <- tkentry(win, width = 50)

  # Create buttons for each of the individual GUIs
  #pd_button <- tkbutton(win, text = "Start PD GUI", command = function() pd(tclvalue(speciesList_entry), tclvalue(speciesFiles_entry)))
  #csp_button <- tkbutton(win, text = "Start CSP GUI", command = function() csp(tclvalue(speciesList_entry), tclvalue(speciesFiles_entry)))
  #cs_button <- tkbutton(win, text = "Start CS GUI", command = function() cs(tclvalue(speciesList_entry), tclvalue(speciesFiles_entry)))

  # Arrange the widgets in a grid
  tkgrid(speciesList_entry, row = 0, column = 0)
  tkgrid(speciesFiles_entry, row = 1, column = 0)
  #tkgrid(pd_button, row = 2, column = 0)
  #tkgrid(csp_button, row = 3, column = 0)
  #tkgrid(cs_button, row = 4, column = 0)
  
  # Set the window title
  tkwm.title(win, "Meta GUI")
}