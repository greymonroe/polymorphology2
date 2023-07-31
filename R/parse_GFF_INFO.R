#' Parse INFO Column from GFF Object
#'
#' This wrapper function runs `parse_INFO` on every value of the INFO column of a GFF object.
#' It takes a GFF object, parses the INFO column, and adds new columns to the object for all the parsed values.
#' GFF is an object returned by read.GFF() and must contain an INFO column.
#'
#' @param GFF A data frame or data.table representing the GFF object. It must contain an INFO column.
#' @return The GFF object with the INFO column parsed into new columns for all values.
#' @examples
#' GFF <- read.GFF("path/to/your.gff")
#' parsed_GFF <- parse_GFF_INFO(GFF)
#' @import data.table
#' @export
parse_GFF_INFO <- function(GFF) {
  # Check if the INFO column exists in the GFF object
  if (!"INFO" %in% names(GFF)) {
    warning("The provided GFF object does not contain an INFO column. Please make sure to provide a valid GFF object with an INFO column.")
    return(GFF)
  }

  # Apply the parse_INFO function to each value of the INFO column
  parsed_data <- lapply(GFF$INFO, parse_INFO)

  # Row-bind the parsed data into a data table
  parsed_data <- data.table::rbindlist(parsed_data, fill = TRUE)

  # Combine the GFF object with the parsed data
  GFF <- cbind(GFF, parsed_data)

  return(GFF)
}
