#' Motif Hunter in Fasta Object
#'
#' This function hunts for a specified DNA sequence motif in a given fasta object.
#' The fasta object must match the case (upper vs lower) of the motif.
#'
#' @param fasta A fasta object as read by read.fasta from the seqinr package in R. Must be a list.
#' @param motif A DNA sequence. Can be upper or lower case but needs to match the case of the fasta object.
#'
#' @return A data.table with locations of the motif within the fasta object.
#'
#' @examples
#' # Ensure fasta object and motif are properly structured before using this function.
#' # motif_hunter(fasta, motif)
#'
#' @export
motif_hunter <- function(fasta, motif){

  # Load necessary package
  if(!requireNamespace("seqinr", quietly = TRUE)){
    stop("Package 'seqinr' is required. Please install it.")
  }

  # Check if fasta object is a list
  if(!is.list(fasta)){
    stop("fasta object must be a list.")
  }

  chrs <- lapply(names(fasta), function(x){
    paste(fasta[[x]], collapse = "")
  })

  allhits <- rbindlist(lapply(names(chrs), function(chr){
    hits <- stringr::str_locate_all(chrs[[chr]], motif)
    hits <- data.table(hits[[1]])
    hits$CHROM <- as.character(chr)
    hits$START <-hits$start
    hits$STOP <-hits$end
    return(hits[,c("CHROM","START","STOP"), with=F])
  }))

  allhits$ID = seq_len(nrow(allhits))
  setkey(allhits, "CHROM","START","STOP")

  return(allhits)
}
