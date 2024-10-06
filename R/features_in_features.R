#' Overlap of Genomic Features
#'
#' This function calculates the overlap between two sets of genomic features.
#' It allows the user to compute counts of overlaps, determine whether any overlap exists,
#' and calculate the mean, sum, length, or sum multiplied by length of a specified value in the second genomic feature set.
#'
#' @param features A `data.table` containing the first set of genomic features.
#'                 It must include the columns "CHROM", "START", "STOP", and "ID".
#' @param features2 A second `data.table` containing the second set of genomic features.
#'                  It must include the columns "CHROM", "START", and "STOP".
#' @param mode The mode of calculation: "counts" for number of overlaps,
#'             "any" for any overlaps (TRUE/FALSE), "mean" for mean of a specified column in the features2 data.table,
#'             "sum" for the sum of a specified column, "length" for the total length of the overlap,
#'             "sumxlength" for the sum of a specified column multiplied by the length of the overlap.
#' @param value The name of a column in the features2 data.table to be used when mode = "mean", "sum", or "sumxlength".
#'
#' @return A `data.table` with the results of the calculation for each feature in the first set.
#'         The result column is named after the 'value' parameter if 'mode' is 'mean', 'sum', or 'sumxlength', or 'overlaps' if 'mode' is 'counts' or 'any'.
#'
#' @examples
#' # Ensure features and features2 data.tables are properly structured before using this function.
#' # Example:
#' # features <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(1, 100, 200), STOP = c(50, 150, 250), ID = c("feat1", "feat2", "feat3"))
#' # features2 <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(25, 125, 225), STOP = c(75, 175, 275))
#' # features_in_features(features, features2, mode = "counts")
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[data.table]{foverlaps}}
#'
#' @importFrom data.table data.table foverlaps setkey
#'
#' @keywords genomics

features_in_features <- function(features, features2, mode, value = NULL) {

  # Check if features and features2 are data.tables
  if (!is.data.table(features) | !is.data.table(features2)) {
    stop("Both 'features' and 'features2' must be data.tables.")
  }

  # Check if necessary columns exist in input data.tables
  required_cols <- c("CHROM", "START", "STOP", "ID")
  if (length(setdiff(required_cols, colnames(features))) > 0 |
      length(setdiff(c("CHROM", "START", "STOP"), colnames(features2))) > 0) {
    stop("The 'features' data.table must have 'CHROM', 'START', 'STOP', and 'ID' columns. 'features2' data.table must have 'CHROM', 'START', and 'STOP' columns.")
  }

  # Check if value is NULL when mode is "sum" or "sumxlength"
  if (mode %in% c("sum", "sumxlength") && is.null(value)) {
    stop("When mode is 'sum' or 'sumxlength', the 'value' parameter must be specified to indicate which column in 'features2' to calculate the sum.")
  }

  # Check if ID column in the features data.table contains unique values
  if (anyDuplicated(features[["ID"]]) > 0) {
    warning("The 'ID' column in the 'features' data.table must contain unique values.")
  }

  # If 'features' contains 'ID', remove it
  if ("ID" %in% colnames(features2)) {
    features2[, ID := NULL]
  }

  # Remove 'LENGTH' column from 'features' if it exists
  if ("LENGTH" %in% names(features)) {
    features$LENGTH <-NULL
  }


  # Add ID column to the features2 data.table
  features2[, feature2_ID := seq_len(.N)]

  # Calculate length in features2
  features2[, LENGTH := STOP - START]

  features$CHROM<-as.character(features$CHROM)
  features2$CHROM<-as.character(features2$CHROM)

  # Check for missing CHROM values in features and features2
  missing_chroms_in_features2 <- setdiff(unique(features$CHROM), unique(features2$CHROM))
  missing_chroms_in_features <- setdiff(unique(features2$CHROM), unique(features$CHROM))

  if (length(missing_chroms_in_features2) > 0) {
    warning(paste("CHROM values in 'features' not found in 'features2':", paste(missing_chroms_in_features2, collapse = ", ")))
  }

  if (length(missing_chroms_in_features) > 0) {
    warning(paste("CHROM values in 'features2' not found in 'features':", paste(missing_chroms_in_features, collapse = ", ")))
  }

  # Set keys for efficient joining
  setkey(features, CHROM, START, STOP)
  setkey(features2, CHROM, START, STOP)

  # Get all unique chromosomes
  unique_chroms <- unique(c(features$CHROM, features2$CHROM))

  # Initiate an empty list to store results for each chromosome
  results_list <- vector("list", length(unique_chroms))

  # Loop over each chromosome
  for (i in seq_along(unique_chroms)) {
    chrom <- unique_chroms[i]
    message(paste("Processing chromosome:", chrom))  # print statement
    features_chrom <- features[CHROM == chrom]
    features2_chrom <- features2[CHROM == chrom]

    # Convert the column to numeric if mode is "sum" or "sumxlength"
    if (mode %in% c("sum", "sumxlength")) {
      features2_chrom[[value]] <- as.numeric(features2_chrom[[value]])
      if (any(is.na(features2_chrom[[value]]))) {
        warning("NAs introduced by coercion to numeric")
      }
    }

    # Perform overlap operation
    overlaps <- foverlaps(features_chrom, features2_chrom)

    # Perform calculations according to the specified mode
    if (mode == "counts") {
      result <- overlaps[, .(counts = sum(!is.na(feature2_ID))), by = .(ID)]
    } else if (mode == "any") {
      result <- overlaps[, .(any = any(!is.na(feature2_ID))), by = .(ID)]
    } else if (mode == "mean" || mode == "sum" || mode == "sumxlength") {
      if (is.null(value)) {
        stop("A 'value' column name must be specified when mode = 'mean', 'sum', or 'sumxlength'.")
      }
      if (!(value %in% colnames(features2_chrom))) {
        stop("The specified 'value' column does not exist in the 'features2' data.table.")
      }
      if (mode == "mean") {
        result <- overlaps[, .(mean = mean(get(value), na.rm = TRUE)), by = .(ID)]
      } else if (mode == "sum") {
        result <- overlaps[, .(sum = sum(get(value), na.rm = TRUE)), by = .(ID)]
      } else if (mode == "sumxlength") {
        result <- overlaps[, .(sumxlength = sum(.SD[[1]], na.rm = TRUE) * sum(LENGTH, na.rm = TRUE)),
                           by = .(ID), .SDcols = value]
      }
    } else if (mode == "length") {
      result <- overlaps[, .(length = sum(LENGTH, na.rm = TRUE)), by = .(ID)]
    } else {
      stop("Invalid 'mode'. Must be one of 'counts', 'any', 'mean', 'sum', 'length', 'sumxlength'.")
    }
    results_list[[i]] <- result
  }

  # Combine results from all chromosomes
  result <- rbindlist(results_list)

  return(result)
}
