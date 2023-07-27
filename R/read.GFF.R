#' Read a GFF Formatted File
#'
#' This function reads a GFF formatted file and assigns the column names.
#' It also checks that the object has the expected number of columns and
#' that START and STOP are numeric.
#'
#' @param file_path The path to the GFF formatted file to read in.
#'
#' @return A data.table with named columns from the GFF file.
#'
#' @export
#' @examples
#' # Ensure the GFF file is in the correct directory before using this function.

read.GFF <- function(file_path) {
  # Check if file exists
  if(!file.exists(file_path)){
    stop("File does not exist")
  }

  # Read the file
  dt <- data.table::fread(file_path, header = FALSE)
  # Check the number of columns
  if(ncol(dt) != 9){
    stop("Unexpected number of columns in the GFF file")
  }

  # Assign column names
  colnames(dt) <- c("CHROM","SOURCE","TYPE","START","STOP","SCORE","DIRECTION","PHASE","INFO")

  # Check if START and STOP columns are numeric
  if(!is.numeric(dt$START) | !is.numeric(dt$STOP)){
    stop("START and STOP columns should be numeric")
  }

  return(dt)
}
