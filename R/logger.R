library(futile.logger)

# Set up the logger
futile.logger::flog.threshold(DEBUG)
futile.logger::flog.appender(appender.file("digestR.log"))

#' Log a Message
#'
#' This function logs a given message using the `futile.logger` package. 
#' Multiple arguments are concatenated into a single log message.
#'
#' @param ... One or more values (of any type). Multiple values are concatenated 
#'   into a single message string.
#'
#' @return NULL. The primary side effect is the logging of the provided message.
#' @importFrom futile.logger flog.info
#'
#' @examples
#' # Log a simple message
#' log_message("This is a test log message.")
#'
#' # Log a message with multiple components
#' log_message("Value of x:", x, "and value of y:", y)
log_message <- function(...){
  # Store all arguments in a list
  args <- list(...)
  
  # Convert all arguments to character and collapse them into a single string
  message <- paste(sapply(args, as.character), collapse = " ")
  
  # Log the message with futile.logger
  futile.logger::flog.info(message)
}
