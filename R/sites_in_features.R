#' Identify Overlapping Sites in Features
#'
#' This function identifies overlapping sites from a provided data.table within given genomic features.
#' The function expects "CHROM", "START", "STOP", and "ID" columns in the features data.table, and
#' "CHROM", "POS" columns in the sites data.table. The 'mode' parameter controls the method
#' of summarizing overlaps: 'counts' gives the total count, 'any' gives a Boolean value
#' indicating the presence of any overlap, and 'mean' calculates the mean of a specified 'value' column.
#'
#' @param features A data.table containing genomic features. Must include "CHROM", "START", "STOP", and "ID" columns.
#' @param sites A data.table of sites (mutations). Must include "CHROM", "POS" columns.
#' @param mode A character string indicating the mode of overlap calculation. Must be one of "counts", "any", or "mean".
#' @param value A character string indicating the column name in 'sites' data.table to calculate the mean. Only required when mode = "mean".
#'
#' @return A data.table summarizing the overlaps for each feature (identified by "ID"), according to the specified mode.
#'
#' @examples
#' # Ensure that features and sites data.tables are properly structured before using this function.
#' # Example:
#' # features <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(1, 100, 200), STOP = c(50, 150, 250), ID = c("feat1", "feat2", "feat3"))
#' # sites <- data.table(CHROM = c("chr1", "chr1", "chr2"), POS = c(25, 125, 225))
#' # sites_in_features(features, sites, mode = "counts")
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[data.table]{foverlaps}}
#'
#' @importFrom data.table data.table foverlaps setkey
#'
#' @keywords genomics

sites_in_features <- function(features, sites, mode, value = NULL) {

  # Check if features and sites are data.tables
  if (!is.data.table(features) | !is.data.table(sites)) {
    stop("Both 'features' and 'sites' must be data.tables.")
  }

  # Check if necessary columns exist in input data.tables
  required_cols <- c("CHROM", "START", "STOP", "ID")
  if (length(setdiff(required_cols, colnames(features))) > 0) {
    stop(paste("The 'features' data.table must have", paste(required_cols, collapse = ", "), "columns."))
  }
  if (length(setdiff(c("CHROM", "POS"), colnames(sites))) > 0) {
    stop("The 'sites' data.table must have 'CHROM' and 'POS' columns.")
  }

  # Check if ID column in the features data.table has unique values
  if (anyDuplicated(features$ID) > 0) {
    stop("The 'ID' column in the 'features' data.table must contain unique values.")
  }

  # Check and remove 'ID' column from 'sites' if exists
  if ("ID" %in% colnames(sites)) {
    sites$ID<-NULL
  }

  # Check for CHROM values not found in the other data.table
  if (length(unmatched <- setdiff(unique(sites$CHROM), unique(features$CHROM))) > 0) {
    warning("The following CHROM values in 'sites' are not found in 'features': ",
            paste(unmatched, collapse = ", "), ".")
  }
  if (length(unmatched <- setdiff(unique(features$CHROM), unique(sites$CHROM))) > 0) {
    warning("The following CHROM values in 'features' are not found in 'sites': ",
            paste(unmatched, collapse = ", "), ".")
  }

  # Add START and STOP columns to the sites data.table
  sites[, c("START", "STOP") := .(POS, POS)]

  features$CHROM<-as.character(features$CHROM)
  sites$CHROM<-as.character(sites$CHROM)

  # Set keys for efficient joining
  setkey(features, CHROM, START, STOP)
  setkey(sites, CHROM, START, STOP)

  # Perform overlap operation
  overlaps <- foverlaps(features, sites)

  # Perform calculations according to the specified mode
  if (mode == "counts") {
    result <- overlaps[, .(counts = sum(!is.na(POS))), by = "ID"]
  } else if (mode == "any") {
    result <- overlaps[, .(any = any(!is.na(POS))), by = "ID"]
  } else if (mode == "mean") {
    if (is.null(value)) {
      stop("A 'value' column name must be specified when mode = 'mean'.")
    }
    if (!(value %in% colnames(sites))) {
      stop("The specified 'value' column does not exist in the 'sites' data.table.")
    }
    result <- overlaps[, .(mean = mean(get(value), na.rm = TRUE)), by = "ID"]
  }

  return(result)
}
