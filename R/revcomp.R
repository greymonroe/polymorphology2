#' Reverse Complement of a DNA Sequence
#'
#' This function returns the reverse complement of a given DNA sequence (motif).
#' It utilizes the `comp` function from the `seqinr` package to generate the complement.
#'
#' @param motif A character string representing the DNA sequence for which the reverse complement is desired.
#'
#' @return A character string representing the reverse complement of the input DNA sequence.
#'
#' @examples
#' revcomp("ATGC")
#' # Returns "GCAT"
#'
#' @importFrom seqinr comp
#'
#' @export
revcomp <- function(motif) {
  toupper(paste0(comp(rev(unlist(strsplit(motif, split = "")))), collapse=""))
}
