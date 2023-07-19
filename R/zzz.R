.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
}

.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("The 'biomaRt' package is required but is not installed. ",
            "Please install it with BiocManager::install('biomaRt').")
  }
}