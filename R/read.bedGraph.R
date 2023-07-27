#' Read bedGraph files
#'
#' This function reads in a bedGraph file, assigns standard column names
#' and checks for the correct format.
#'
#' @param file A character string specifying the path to the bedGraph file.
#'
#' @return A `data.table` with the bedGraph data, where columns are named "CHROM", "START", "STOP", and "DEPTH".
#'
#' @examples
#' # Assuming 'example.bedGraph' is a correctly formatted bedGraph file in your working directory:
#' # data <- read.bedGraph('example.bedGraph')
#'
#' @export
#'
#' @importFrom data.table fread
#'
#' @keywords fileIO genomics

read.bedGraph <- function(file) {
  # Read the file
  data <- fread(file)

  # Check that the file has exactly four columns
  if (ncol(data) != 4) {
    stop("The input file must have exactly four columns.")
  }

  # Rename the columns
  colnames(data) <- c("CHROM", "START", "STOP", "DEPTH")

  # Check that START, STOP, and DEPTH are numeric
  if (!is.numeric(data$START) | !is.numeric(data$STOP) | !is.numeric(data$DEPTH)) {
    stop("The 'START', 'STOP', and 'DEPTH' columns must be numeric.")
  }

  return(data)
}
