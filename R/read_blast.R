#' Read BLAST Output
#'
#' This function reads BLAST output in format 6 and returns the results in a data.table format with appropriate column names.
#'
#' @param file A string. The file path corresponding to the output of a BLAST search in fmt 6.
#'
#' @return A data.table object containing BLAST results with appropriate column names.
#'
#' @examples
#' # Assuming 'sample_blast_output.txt' is a BLAST output in fmt 6
#' blast_results <- read_blast('sample_blast_output.txt')
#'
#' @import data.table
#' @export
read_blast <- function(file) {
  # Ensure data.table library is loaded
  require(data.table)

  # Read the file
  blast <- fread(file, header = FALSE)

  # Check if blast is empty, and if so, create a placeholder table with NA values
  if (nrow(blast) == 0) {
    blast <- data.table(matrix(NA, ncol = 12, nrow = 1))
  }

  # Assign column names
  colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore")

  return(blast)
}
