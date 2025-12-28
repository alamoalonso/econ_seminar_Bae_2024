#' Logging Utilities
#'
#' Simple logging system for debugging and tracing.
#'
#' @name utils_logging
NULL

#' Log an informational message
#'
#' @param msg Character string to log
#' @param config List containing debug flag
#' @export
log_info <- function(msg, config = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[INFO %s] %s\n", timestamp, msg))
}

#' Log a warning message
#'
#' @param msg Character string to log
#' @param config List containing debug flag
#' @export
log_warn <- function(msg, config = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[WARN %s] %s\n", timestamp, msg))
}

#' Log a debug message (only if debug flag is TRUE)
#'
#' @param msg Character string to log
#' @param config List containing debug flag (config$debug)
#' @export
log_debug <- function(msg, config = NULL) {
  if (!is.null(config) && isTRUE(config$debug)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[DEBUG %s] %s\n", timestamp, msg))
  }
}

#' Check if current time index should be traced
#'
#' @param t_idx Current time index
#' @param config List containing trace_origins vector
#' @return Logical indicating whether to trace
#' @export
should_trace <- function(t_idx, config = NULL) {
  if (is.null(config) || is.null(config$trace_origins)) {
    return(FALSE)
  }
  t_idx %in% config$trace_origins
}
