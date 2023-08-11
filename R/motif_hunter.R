#' Motif Hunter in Fasta Object
#'
#' This function hunts for a specified DNA sequence motif in a given fasta object.
#' The fasta object must match the case (upper vs lower) of the motif.
#'
#' When ambiguous nucleotides (IUPAC codes) are included in the motif, the function
#' expands the motif into all possible combinations based on these codes and searches
#' for each combination in the fasta object. If the `reverse_complement` parameter is set to TRUE,
#' the function will also search for the reverse complement of each motif. The result
#' will include an additional column, 'MOTIF', to specify which version of the expanded
#' motif (or its reverse complement) was found.
#'
#' @param fasta A fasta object as read by read.fasta from the seqinr package in R. Must be a list.
#' @param motif A DNA sequence, which can include IUPAC ambiguous nucleotide codes. Automatically converted to upper case.
#' @param reverse_complement Logical. If TRUE, the function will also search for the reverse complement of the motif.
#'
#' @return A data.table with locations of the motif within the fasta object. If ambiguous nucleotides are used in the motif, an additional column 'MOTIF' will specify which version of the motif was found.
#'
#' @examples
#' # Ensure fasta object and motif are properly structured before using this function.
#' # motif_hunter(fasta, motif, reverse_complement = TRUE)
#'
#' @export
#'
#' @seealso \code{\link[data.table]{data.table}}, \code{\link[seqinr]{read.fasta}}, IUPAC nucleotide code
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom stringr str_locate_all
#'
#' @keywords genomics
motif_hunter <- function(fasta, motif, reverse_complement){

  # Load necessary package
  if(!requireNamespace("seqinr", quietly = TRUE)){
    stop("Package 'seqinr' is required. Please install it.")
  }

  # Check if fasta object is a list
  if(!is.list(fasta)){
    stop("fasta object must be a list.")
  }

  # Ambiguous nucleotide codes
  ambiguous_codes <- list(
    R = c("A", "G"),
    Y = c("C", "T"),
    S = c("G", "C"),
    W = c("A", "T"),
    K = c("G", "T"),
    M = c("A", "C"),
    B = c("C", "G", "T"),
    D = c("A", "G", "T"),
    H = c("A", "C", "T"),
    V = c("A", "C", "G"),
    N = c("A", "C", "G", "T")
  )

  # Expand motif into all possible combinations based on ambiguous nucleotides
  expand_motif <- function(motif, ambiguous_codes) {
    if (length(grep("[RYSWKMBDHVN]", motif)) == 0) {
      return(list(motif))
    } else {
      first_ambig <- regmatches(motif, regexpr("[RYSWKMBDHVN]", motif))[[1]]
      motifs <- lapply(ambiguous_codes[[first_ambig]], function(nucleotide) {
        return(gsub(first_ambig, nucleotide, motif, fixed = TRUE))
      })
      return(unlist(lapply(motifs, function(m) expand_motif(m, ambiguous_codes))))
    }
  }

  motifs <- expand_motif(toupper(motif), ambiguous_codes)

  # If reverse_complement is TRUE, add reverse complements to the motifs list
  if (reverse_complement) {
    motifs <- c(motifs, lapply(motifs, revcomp))
  }

  chrs <- lapply(names(fasta), function(x){
    toupper(paste(fasta[[x]], collapse = ""))
  })
  names(chrs) <- names(fasta)

  allhits <- rbindlist(lapply(motifs, function(motif) {
    rbindlist(lapply(names(chrs), function(chr){
      hits <- stringr::str_locate_all(chrs[[chr]], motif)
      hits <- data.table(hits[[1]])
      hits$CHROM <- as.character(chr)
      hits$START <- hits$start
      hits$STOP <- hits$end
      hits$MOTIF <- motif
      return(hits[, c("CHROM", "START", "STOP", "MOTIF"), with = F])
    }))
  }))

  allhits$ID = seq_len(nrow(allhits))
  setkey(allhits, "CHROM", "START", "STOP")

  return(allhits)
}
