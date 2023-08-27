#' Calculate ChIP Enrichment for Genomic Features
#'
#' This function calculates the ChIP enrichment for a set of genomic features using ChIP and input bedGraph files.
#' It reads the bedGraph files, calculates the overlaps with the features, and calculates the enrichment.
#'
#' @param features A `data.table` containing genomic features. It must include the columns "CHROM", "START", "STOP", and "ID".
#' @param chipfile The name of the ChIP bedGraph file.
#' @param inputfile The name of the input bedGraph file.
#'
#' @return A `data.table` with the calculated enrichment for each feature.
#'
#' @examples
#' # Ensure features data.table is properly structured before using this function.
#' # Example:
#' # features <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(1, 100, 200), STOP = c(50, 150, 250), ID = c("feat1", "feat2", "feat3"))
#' # features_chip_enrich(features = features, chipfile = "chip.bedGraph", inputfile = "input.bedGraph")
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[data.table]{read.bedGraph}}, \code{\link[data.table]{bedGraph_total}}, \code{\link[data.table]{features_in_features}}
#'
#' @keywords genomics
#'

features_chip_enrich <- function(features, chipfile, inputfile) {

  message("Reading input bedGraph file...")
  input <- read.bedGraph(inputfile)

  message("Reading ChIP bedGraph file...")
  chip <- read.bedGraph(chipfile)

  message("Checking if chromosome names match in all objects...")
  if (!all(unique(features$CHROM) %in% unique(chip$CHROM)) | !all(unique(features$CHROM) %in% unique(input$CHROM))) {
    warning("Chromosome names in 'features', 'chip', and 'input' have mismatches.")
  }

  message("Calculating overlap of features with input...")
  features_input <- features_in_features(features = features, features2 = input, mode = "sumxlength", value="DEPTH")
  colnames(features_input)[2] <- "input"

  message("Calculating overlap of features with ChIP...")
  features_chip <- features_in_features(features = features, features2 = chip, mode = "sumxlength", value="DEPTH")
  colnames(features_chip)[2] <- "chip"

  message("Merging overlap results...")
  merge <- merge(features_chip, features_input)

  message("Calculating total depth of ChIP...")
  chip_total <- bedGraph_total(chip)

  message("Calculating total depth of input...")
  input_total <- bedGraph_total(input)

  message("Calculating enrichment...")
  merge$enrich <- log2((1 + merge$chip) / chip_total) - log2((1 + merge$input) / input_total)

  message("Enrichment calculation completed successfully.")

  return(merge)
}

