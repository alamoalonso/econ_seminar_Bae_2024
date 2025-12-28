#' Data I/O Functions
#'
#' Functions for loading datasets from various sources.
#'
#' @name data_io
NULL

#' Load dataset based on configuration
#'
#' @param config Configuration list containing dataset_id and data_file
#' @return A list with components:
#'   - raw_data: data.frame with date column and predictor columns
#'   - dates: vector of Date objects
#'   - dataset_id: character string identifying the dataset
#' @export
#' @examples
#' config <- config_us_default()
#' dataset <- load_dataset(config)
load_dataset <- function(config) {
  log_info(sprintf("Loading dataset: %s", config$dataset_id), config)

  if (config$dataset_id == "US_FRED") {
    data <- load_fred_md(config)
  } else if (config$dataset_id == "EU_EA_MD_QD") {
    data <- load_ea_md_qd(config)
  } else {
    stop(sprintf("Unknown dataset_id: %s", config$dataset_id))
  }

  log_info(sprintf("Loaded %d observations x %d variables", nrow(data$raw_data), ncol(data$raw_data) - 1), config)

  data
}

#' Load FRED-MD dataset
#'
#' @param config Configuration list
#' @return List with raw_data, dates, dataset_id
#' @importFrom fbi fredmd
load_fred_md <- function(config) {
  if (!file.exists(config$data_file)) {
    stop(sprintf("Data file not found: %s", config$data_file))
  }

  # Load and transform using fbi::fredmd
  fred_transformed <- fbi::fredmd(
    file = config$data_file,
    transform = config$apply_transforms
  )

  # Coerce to plain data.frame and strip 'fredmd' class
  fred_df <- as.data.frame(fred_transformed)
  class(fred_df) <- "data.frame"

  # Apply sample period filter
  fred_df <- fred_df[fred_df$date >= config$sample_start & fred_df$date <= config$sample_end, ]

  if (nrow(fred_df) == 0) {
    stop("No observations remain after filtering by sample_start/sample_end")
  }

  list(
    raw_data = fred_df,
    dates = fred_df$date,
    dataset_id = config$dataset_id
  )
}

#' Load Euro Area EA-MD-QD dataset (author-processed .xlsx output)
#'
#' @param config Configuration list
#' @return List with raw_data, dates, dataset_id
#' @importFrom readxl read_xlsx
load_ea_md_qd <- function(config) {
  if (!file.exists(config$data_file)) {
    stop(sprintf("EA-MD-QD data file not found: %s", config$data_file))
  }
  
  ext <- tolower(tools::file_ext(config$data_file))
  if (ext != "xlsx") {
    stop(sprintf("EA-MD-QD loader expects .xlsx produced by routine_data.py, got: .%s", ext))
  }
  
  df <- readxl::read_xlsx(config$data_file)
  
  # Author output uses 'Date' as index name; accept both
  if ("Date" %in% names(df) && !("date" %in% names(df))) {
    names(df)[names(df) == "Date"] <- "date"
  }
  if (!("date" %in% names(df))) {
    stop("EA-MD-QD file must contain a 'Date' or 'date' column.")
  }
  
  # Robust date parsing
  if (inherits(df$date, c("POSIXct", "POSIXt"))) {
    df$date <- as.Date(df$date)
  } else if (inherits(df$date, "Date")) {
    # keep
  } else if (is.numeric(df$date)) {
    # Excel serial date fallback (1899-12-30 origin is standard for Windows Excel)
    df$date <- as.Date(df$date, origin = "1899-12-30")
  } else {
    # character or other -> try standard parsing
    df$date <- as.Date(df$date)
  }
  
  if (anyNA(df$date)) {
    stop("Date parsing produced NA values. Inspect the 'Date' column in the .xlsx.")
  }
  
  # Ensure date is first column (preprocess_dataset assumes this)
  df <- df[, c("date", setdiff(names(df), "date"))]
  
  # Apply sample period filter
  df <- df[df$date >= config$sample_start & df$date <= config$sample_end, , drop = FALSE]
  if (nrow(df) == 0) stop("No observations remain after filtering by sample_start/sample_end")
  
  list(
    raw_data   = as.data.frame(df),
    dates      = df$date,
    dataset_id = config$dataset_id
  )
}
