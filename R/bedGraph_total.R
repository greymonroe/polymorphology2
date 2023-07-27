#' Calculate Total Depth of a BedGraph
#'
#' This function calculates the total depth of a bedGraph.
#' The total depth is calculated by summing the product of DEPTH and the length of the interval (STOP - START) for each row.
#'
#' @param bed A data.table returned by the function read.bedGraph().
#'            It must include the columns "CHROM", "START", "STOP", and "DEPTH".
#'
#' @return The total depth of the bedGraph.
#'
#' @examples
#' # Ensure bedGraph data.table is properly structured before using this function.
#' # Example:
#' # bed <- data.table(CHROM = c("chr1", "chr1"), START = c(1, 100), STOP = c(50, 150), DEPTH = c(2, 3))
#' # bedGraph_total(bed)
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[data.table]{read.bedGraph}}
#'
#' @keywords genomics
#'
#' @importFrom data.table data.table

bedGraph_total <- function(bed){

  # Check if bed is a data.table
  if (!is.data.table(bed)) {
    stop("'bed' must be a data.table.")
  }

  # Check if necessary columns exist in input data.table
  required_cols <- c("CHROM", "START", "STOP", "DEPTH")
  if (any(setdiff(required_cols, colnames(bed)))) {
    stop("The 'bed' data.table must have 'CHROM', 'START', 'STOP', and 'DEPTH' columns.")
  }

  # Calculate total depth
  total <- sum(bed$DEPTH * (bed$STOP - bed$START))

  return(total)
}
