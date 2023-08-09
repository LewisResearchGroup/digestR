#.onAttach <- function(libname, pkgname) {
#  if (!requireNamespace("biomaRt", quietly = TRUE)) {
#    if (!requireNamespace("BiocManager", quietly = TRUE)) {
#      install.packages("BiocManager")
#    }
#    BiocManager::install("biomaRt", ask = FALSE)
#    message("`biomaRt` was missing and has been installed.")
#  }
#}
