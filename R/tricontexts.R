#' Trinucleotide context
#'
#' This function obtains the trinucleotide context for a given set of mutations.
#' The context is determined by examining the base pairs immediately flanking
#' the mutation site.
#'
#' @param sites A data.table of mutations. It should contain "CHROM", "POS",
#' "REF", and "ALT" columns, representing chromosome, position, reference allele,
#' and alternate allele, respectively.
#' @param fasta A named list representing the genome sequence. Names of the list
#' should correspond to chromosome names, and each item in the list should be a
#' character vector representing the sequence of the corresponding chromosome.
#'
#' @return A vector of contexts for each mutation.
#' @export

tricontexts <- function(sites, fasta) {
  SBS <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

  # Check if all chromosomes in the sites are present in the names of fasta
  chroms_not_found <- setdiff(sites[, CHROM], names(fasta))
  if (length(chroms_not_found) > 0) {
    warning("The following chromosomes from 'sites' are not found in the 'fasta': ",
            paste(chroms_not_found, collapse = ", "))
  }

  # create progress bar
  pb <- txtProgressBar(min = 0, max = nrow(sites), style = 3)

  contexts <- sapply(seq_len(nrow(sites)), function(i) {
    chr <- sites[i, CHROM]
    pos <- as.numeric(sites[i, POS])
    ref <- sites[i, REF]
    alt <- sites[i, ALT]

    if(nchar(alt) != 1 | nchar(ref) != 1){
      return(NA)
    }

    if(!(chr %in% names(fasta))){
      return(NA)
    }

    if(pos>length(fasta[[chr]])){
      stop(paste("row",i,"POS is larger than chromosome length in fasta"))
    }
    up <- toupper(paste0(fasta[[chr]][(pos-1)], collapse = ""))
    down <- toupper(paste0(fasta[[chr]][(pos+1)], collapse = ""))
    context <- paste0(up,"[",ref,">",alt,"]", down)
    mutation <- paste0(up,ref,down,">", up,alt,down)

    if(!context %in% SBS){
      ref <- toupper(comp(ref))
      alt <- toupper(comp(alt))
      up <- toupper(comp(up))
      down <- toupper(comp(down))
      context <- paste0(down,"[",ref,">",alt,"]", up)
      mutation <- paste0(up,ref,down,">", up,alt,down)
    }

    # update progress bar
    setTxtProgressBar(pb, i)

    return(mutation)
  })

  # close progress bar
  close(pb)

  return(contexts)
}
