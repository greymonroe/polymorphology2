#' Parse INFO Field from GFF Format
#'
#' This function takes a single INFO field from the GFF (General Feature Format) file
#' and parses it into a data table. The INFO field consists of a semicolon-separated
#' string of key-value pairs, where each key and value are separated by an equal sign.
#' Note: This is a subfunction run by parse_GFF_INFO.
#'
#' @param INFO A character string representing the INFO field from a GFF file. It can be obtained using read.GFF()$INFO.
#' @return A data.table containing the parsed key-value pairs from the INFO field.
#' @examples
#' info_string <- "key1=value1;key2=value2;key3=value3"
#' parsed_info <- parse_INFO(info_string)
#' @import data.table
#' @export
parse_INFO <- function(INFO) {
  # Split the INFO string by semicolons to separate the key-value pairs
  key_values <- unlist(strsplit(INFO, split = ";"))

  # Extract the values from the key-value pairs
  parsed_item <- sapply(key_values, function(x) unlist(strsplit(x, split = "="))[2])

  # Extract the keys from the key-value pairs
  col_names <- sapply(key_values, function(x) unlist(strsplit(x, split = "="))[1])

  # Create a data table with the parsed key-value pairs
  result <- data.table::data.table(t(parsed_item))

  # Set the column names of the data table to the extracted keys
  data.table::setnames(result, col_names)

  return(result)
}
