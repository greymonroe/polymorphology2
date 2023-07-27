#' Identify Overlapping Features in Sites
#'
#' This function identifies whether sites from a provided data.table overlap with given genomic features.
#' The function expects "CHROM", "START", and "STOP" columns in the features data.table, and
#' "CHROM", "POSITION", and "ID" columns in the sites data.table.
#'
#' @param features A data.table containing genomic features. Must include "CHROM", "START", and "STOP" columns.
#'                 If an "ID" column is present, it will be removed.
#' @param sites A data.table of sites (mutations). Must include "CHROM", "POSITION", and "ID" columns.
#'
#' @return A data.table indicating whether each site (identified by "ID") overlaps with any feature.
#'
#' @examples
#' # Ensure that features and sites data.tables are properly structured before using this function.
#' # Example:
#' # features <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(1, 100, 200), STOP = c(50, 150, 250))
#' # sites <- data.table(CHROM = c("chr1", "chr1", "chr2"), POSITION = c(25, 125, 225), ID = c("site1", "site2", "site3"))
#' # features_in_sites(features, sites)
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[data.table]{foverlaps}}
#'
#' @importFrom data.table data.table foverlaps setkey
#'
#' @keywords genomics

features_in_sites <- function(features, sites) {

  # Check if necessary columns exist in input data.tables
  if (any(setdiff(c("CHROM", "START", "STOP"), colnames(features)))){
    stop("features data.table needs to have CHROM, START, and STOP columns.")
  }
  if (any(setdiff(c("CHROM", "POSITION", "ID"), colnames(sites)))){
    stop("sites data.table needs to have CHROM, POSITION, and ID columns.")
  }

  # If 'features' contains 'ID', remove it
  if ("ID" %in% colnames(features)) {
    warning("'features' data.table contains 'ID'. This column will be removed.")
    features[, ID := NULL]
  }

  # Add feature_ID to features data.table
  features[, feature_ID := seq_len(.N)]

  # Create START and STOP columns in the sites data.table
  sites[, c("START", "STOP") := .(POSITION, POSITION)]

  # Convert CHROM columns to character
  features[, CHROM := as.character(CHROM)]
  sites[, CHROM := as.character(CHROM)]

  # Check for unmatched CHROM values
  unmatched_features <- setdiff(unique(features$CHROM), unique(sites$CHROM))
  unmatched_sites <- setdiff(unique(sites$CHROM), unique(features$CHROM))
  if (length(unmatched_features) > 0) {
    warning(paste("The following CHROM values in 'features' are not found in 'sites':",
                  paste(unmatched_features, collapse=", ")))
  }
  if (length(unmatched_sites) > 0) {
    warning(paste("The following CHROM values in 'sites' are not found in 'features':",
                  paste(unmatched_sites, collapse=", ")))
  }

  # Set keys for efficient joining
  setkey(features, CHROM, START, STOP)
  setkey(sites, CHROM, START, STOP)

  # Identify overlaps
  overlaps <- foverlaps(sites, features)

  # Summarize overlaps for each site
  out <- overlaps[, .(overlaps = any(!is.na(feature_ID))), by = .(ID)]

  return(out)
}
