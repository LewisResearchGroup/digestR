# This script provides a global logging function

# Define the logger function
log_message <- function(...){
  # Store all arguments in a list
  args <- list(...)
  
  # Use do.call to apply the print function to the list of arguments
  do.call(print, args)
  
  # Define the log file
  log_file <- "digestR.log"
  
  # Convert all arguments to character and collapse them into a single string
  message <- paste(sapply(args, as.character), collapse = " ")
  
  # Write the log message to the file
  write(paste(Sys.time(), ": ", message, "\n", sep = ""),
        file = log_file, 
        append = TRUE)
}