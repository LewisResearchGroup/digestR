#' @importFrom grDevices X11 cm.colors col2rgb contourLines dev.cur dev.list dev.new dev.next dev.off dev.set rgb topo.colors
NULL

#' @importFrom graphics .filled.contour abline axis box contour grconvertX grconvertY image legend lines locator par persp points polygon rect segments text title
NULL

#' @importFrom stats approx fivenum median na.exclude na.omit sd smooth.spline
NULL

#' @importFrom utils Rprof available.packages browseURL capture.output compareVersion contrib.url download.file flush.console install.packages installed.packages packageDescription read.csv read.table select.list tail tar untar update.packages write.csv write.table
NULL


maybe_create_directory <- function(dir) {
  # Check if the directory exists
  if (!dir.exists(dir)) {
    # If the directory does not exist, create it
    dir.create(dir, recursive = TRUE)  # The recursive = TRUE argument means that it will create any necessary parent directories as well
    log_message('Created directory:', dir)
  }
}


set_directory <- function(dir) {
  maybe_create_directory(dir)
  # Set the working directory
  setwd(dir)
  log_message('Working directory:', getwd())
}


recursive_delete <- function(dir) {
  unlink(dir, recursive = TRUE, force = TRUE)
}