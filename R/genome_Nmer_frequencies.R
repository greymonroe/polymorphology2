
#' Calculate N-mer Frequencies in Genomic Sequences
#'
#' This function calculates the frequencies of N-mers (specified by the user) in genomic sequences
#' provided as a FASTA object. It is designed to efficiently process large genomic sequences by
#' optionally limiting the analysis to the first `stop` positions of each chromosome. The results
#' are summarized in a data.table showing the counts of each N-mer across all processed sequences.
#'
#' @param genome A FASTA object as read by `seqinr::read.fasta()`. This should be a list where
#'        each element represents a chromosome or sequence, with names corresponding to chromosome
#'        or sequence identifiers.
#' @param Nmer An integer specifying the size of the N-mer to analyze.
#' @param stop An optional integer. If provided, the analysis is limited to the first `stop`
#'        positions of each sequence. By default (`NULL`), the entire length of each sequence is
#'        considered.
#'
#' @return A data.table with two columns: `context`, representing the unique N-mers identified across
#'         all sequences, and `N`, representing the sum of occurrences of each N-mer. This table
#'         provides a comprehensive overview of N-mer frequencies within the specified genomic regions.
#'
#' @examples
#' # Assuming `fasta_object` is read using seqinr::read.fasta()
#' Nmer_counts <- calculate_Nmer_frequencies(fasta_object, Nmer = 3, stop = 10000)
#'
#' @import data.table
#' @export
genome_Nmer_frequencies <- function(genome, Nmer, stop = NULL) {

  #stop <- if(is.null(stop)) length(genome[[n]]) else min(stop, length(genome[[n]]))

  genome_features <- rbindlist(lapply(names(genome), function(n) {
    chr <- genome[[n]]
    out <- data.table(CHROM = n, START = 1, STOP = if(is.null(stop)) length(chr) else stop)
    return(out)
  }))

  genome_features$ID <- 1:nrow(genome_features)
  genome_features_trimer <- Nmerfrequency(genome_features, genome, Nmer = Nmer, mode = "counts")

  # Assuming Nmerfrequency returns a data.table with N-mer counts, aggregate to get total counts per N-mer
  genome_features_trimer_total <- rbindlist(lapply(1:length(genome_features_trimer), function(x){
    data.table(context=names(genome_features_trimer)[x], N=sum(genome_features_trimer[,..x]))
  }))
  genome_features_trimer_total<-genome_features_trimer_total[!grepl("ID|N|V|Y|W|R|M|K|S", context)]
  genome_features_trimer_total$TRI<-genome_features_trimer_totalgenom_tri$context
  genome_features_trimer_total<-genome_features_trimer_total[1:64]
  return(genome_features_trimer_total)
}
