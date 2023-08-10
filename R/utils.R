#' @importFrom grDevices X11 cm.colors col2rgb contourLines dev.cur dev.list dev.new dev.next dev.off dev.set rgb topo.colors
NULL

#' @importFrom graphics .filled.contour abline axis box contour grconvertX grconvertY image legend lines locator par persp points polygon rect segments text title
NULL

#' @importFrom stats approx fivenum median na.exclude na.omit sd smooth.spline
NULL

#' @importFrom utils Rprof available.packages browseURL capture.output choose.dir compareVersion contrib.url download.file flush.console install.packages installed.packages packageDescription read.csv read.table select.list tail tar untar update.packages winMenuAdd winMenuAddItem winMenuNames write.csv write.table
NULL


maybeCreateDirectory <- function(dir) {
  # Check if the directory exists
  if (!dir.exists(dir)) {
    # If the directory does not exist, create it
    dir.create(dir, recursive = TRUE)  # The recursive = TRUE argument means that it will create any necessary parent directories as well
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