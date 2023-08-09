library(futile.logger)

# Set up the logger
futile.logger::flog.threshold(DEBUG)
futile.logger::flog.appender(appender.file("digestR.log"))

# Define the logger function
log_message <- function(...){
  # Store all arguments in a list
  args <- list(...)
  
  # Convert all arguments to character and collapse them into a single string
  message <- paste(sapply(args, as.character), collapse = " ")
  
  # Log the message with futile.logger
  futile.logger::flog.info(message)
}
