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
  # Split the INFO string by semicolons to get potential key-value pairs
  INFO<-gsub("\\(source=","source", INFO)
  potential_key_values <- unlist(strsplit(INFO, split = ";"))

  # Initialize variables for keys and values
  keys <- vector("character")
  values <- vector("character")

  # Iterate through potential key-value pairs
  for (item in potential_key_values) {
    if (grepl("=", item)) {
      # Split by equal sign to extract key and value
      kv <- unlist(strsplit(item, split = "="))
      keys <- c(keys, kv[1])
      values <- c(values, kv[2])
    } else {
      # Append item to the last value if no equal sign is found
      values[length(values)] <- paste(values[length(values)], item, sep = ";")
    }
  }

  # Create a data table with the parsed keys and values
  result <- data.table::data.table(t(values))

  # Set the column names of the data table to the extracted keys
  data.table::setnames(result, keys)

  return(result)
}
