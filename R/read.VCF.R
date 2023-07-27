
#' Read and parse a VCF file
#'
#' This function uses the read.vcfR function from the vcfR package to read and parse a VCF file.
#' It allows for input from different callers including Strelka2, HaplotypeCaller, and PBSV.
#'
#' @param file Path to the VCF file to be read.
#' @param caller Name of the caller used. Options are "strelka2", "haplotypecaller", "pbsv","other".
#'
#' @return A data.table containing the fixed fields and genotype information from the VCF file.
#' If the caller is "strelka2", an additional EVS column is added to the output.
#'
#' @examples
#' # Ensure vcfR package is installed before using this function.
#' # read.VCF("path_to_your_file.vcf", "strelka2")
#'
#' @export
read.VCF <- function(file, caller){

  # Ensure vcfR package is installed
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package 'vcfR' needed for this function to work. Please install it.")
  }

  # Read the VCF file
  dat <- vcfR::read.vcfR(file, verbose = F)

  # Create a data.table with fixed fields and genotype information
  variants <- cbind(data.table(dat@fix), data.table(dat@gt))

  # If the caller is Strelka2, add the EVS information to the output
  if(caller == "strelka2"){
    variants$EVS <- as.numeric(gsub(".+SomaticEVS=", "", variants$INFO))
  }

  variants$POS<-as.numeric(variants$POS)

  return(variants)
}
