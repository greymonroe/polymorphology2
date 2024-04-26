#' Annotate Variants Near Homopolymers
#'
#' This function identifies and annotates variants that are in or near homopolymer runs in the genome.
#' Variants are annotated based on whether they are located within a certain distance (specified by 'dist')
#' from a homopolymer run of a certain size (specified by 'size').
#'
#' @param sites A data.table with CHROM, POS, REF, ALT columns representing the genomic sites.
#' @param homopolymers A data.table returned from the `make_homopolymer()` function, representing homopolymer regions by CHROM.
#' @param size The minimum length of homopolymers to consider.
#' @param dist The maximum distance a variant can be from a homopolymer to be considered a neighbor.
#'
#' @return A data.table with a new column, 'homopolymer_neighbor', which indicates whether each site is a neighbor to a homopolymer.
#'
#' @examples
#' # Assume 'sites' and 'homopolymers' are pre-defined
#' # sites <- data.table( ... )
#' # homopolymers <- make_homopolymer( ... )
#' # Annotate sites
#' sites <- homopolymer_var_annotate(sites, homopolymers, 5, 10)
#'
#' @export
#'
#' @seealso \code{\link[data.table]{foverlaps}}, \code{\link[data.table]{rbindlist}}, \code{\link[data.table]{setkey}}
#'
#' @keywords genomics


homopolymer_var_annotate <- function(sites, homopolymers, size, dist){

  # Ensure that 'sites' has the correct columns
  if (!all(c("CHROM", "POS", "REF", "ALT") %in% names(sites))) {
    stop("'sites' must have 'CHROM', 'POS', 'REF', and 'ALT' columns.")
  }

  # Ensure that 'homopolymers' has the correct structure
  if (!all(c("CHROM", "START", "STOP", "LENGTH","BP" ) %in% names(homopolymers))) {
    stop("'homopolymers' must have 'CHROM', 'START', 'STOP',  'LENGTH' and 'BP' columns.")
  }


  # Add columns to homopolymers for start minus 'dist' and end plus 'dist'
  homopolymers$startminus <- homopolymers$START - dist
  homopolymers$endplus <- homopolymers$STOP + dist
  homopolymers$CHROM<-as.character(homopolymers$CHROM)


  # Set keys for overlap matching
  data.table::setkey(homopolymers, CHROM, startminus, endplus)
  sites$START<-sites$POS
  sites$STOP<-sites$POS
  sites$CHROM<-as.character(sites$CHROM)
  data.table::setkey(sites, CHROM, START, STOP)

  # Identify overlaps between sites and homopolymers

  overlap_sums<-rbindlist(lapply(unique(sites$CHROM), function(c){
    message("CHROM", c)
    overlap <- data.table::foverlaps(sites[CHROM==c], homopolymers[CHROM==c])

    # Summarize overlap information
    overlap_sum <- overlap[toupper(ALT) == toupper(BP), .(homopolymer_neighbor = sum(!is.na(BP)) > 0), by = ID]

  }))


  # Annotate sites with overlap information
  sites$homopolymer_neighbor = overlap_sums$homopolymer_neighbor[match(sites$ID, overlap_sums$ID)]
  sites$homopolymer_neighbor[is.na(sites$homopolymer_neighbor)]<-0

  return(sites[,.(ID, HP=homopolymer_neighbor)])
}
