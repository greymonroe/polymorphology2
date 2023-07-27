#' Identify Homopolymer Runs
#'
#' This function is a wrapper for the `homopolymer_read` function. It is used
#' to identify runs of homopolymers in a genome sequence.
#'
#' @param fasta A named list representing the genome sequence.
#' @param size The minimum length of homopolymers to identify. This parameter
#'             is passed to the `homopolymer_read` function.
#'
#' @return A list of identified homopolymers for each sequence in the genome.
#'
#' @examples
#' # Define a list representing a genome sequence
#' genome_seq <- list("seq1" = c("A", "A", "A", "T", "G", "G", "G", "G"),
#'                    "seq2" = c("C", "C", "C", "A", "A", "T", "T", "T"))
#' # Call the function to identify homopolymers
#' homopolymers <- make_homopolymer(genome_seq, 2)
#'
#' @export
#'
#' @seealso \code{\link[homopolymer_read]{homopolymer_read}}
#'
#' @keywords genomics

make_homopolymer <- function(fasta, size) {
  # Split each sequence into individual nucleotides
  seqs<-lapply(fasta, function(x) unlist(x[1:length(x)]))
  homopolymers<-lapply(seqs, function(x) homopolymer_read(x, size))

  # Ensure homopolymers meet the specified size requirement and add CHROM column
  homopolymers <- data.table::rbindlist(lapply(1:length(homopolymers), function(x) {
    out <- homopolymers[[x]][LENGTH >= size]
    out$CHROM <- names(homopolymers)[x]
    return(out)
  }))

  homopolymers<-homopolymers[,.(CHROM, START, STOP, LENGTH, BP)]
  return(homopolymers)
}
