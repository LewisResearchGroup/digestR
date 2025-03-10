if (requireNamespace("digestR", quietly = TRUE)) {
  cat("digestR is already installed.\n")
  return()  # Exit the script if the package is installed
}
​
sink("output.log", append = TRUE, type = "output")
​
#Sys.sleep(1)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = Sys.getenv('R_LIBS_USER'), repos = 'https://cloud.r-project.org', dependencies = TRUE, INSTALL_opts = '--no-lock')
​
​
​
​
cat("checkpoint 1: BiocManager\n")
​
​
​
if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("biomaRt", quietly = TRUE)) {
  # If either package is not installed, install them
  BiocManager::install(c("AnnotationDbi", "biomaRt"))
} else {
  cat("Both AnnotationDbi and biomaRt are already installed.\n")
}
​
cat("checkpoint 2: AnnotationDbi and biomaRt\n")
#Sys.sleep(1)
​
​
# Check if devtools is installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  # If devtools is not installed, install it from CRAN
  install.packages("devtools", repos = "https://cran.rstudio.com/")
} else {
  cat("devtools is already installed.\n")
​
}
​
#Sys.sleep(1)
​
library(devtools)
​
cat("checkpoint 3: devtools\n")
#Sys.sleep(1)
​
​
if (!requireNamespace("digestR", quietly = TRUE)) {
  devtools::install_github("LewisResearchGroup/digestR")
} else {
  cat("digestR is installed.\n")
}
​
​
cat("checkpoint 4: digestR\n")
Sys.sleep(1)
​
sink()
cat("installation successful.\n")
​
#library(digestR)
