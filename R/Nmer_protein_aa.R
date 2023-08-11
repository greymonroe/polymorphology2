#' Nmer_protein_aa
#'
#' This function calculates the frequency or proportion of N-mer amino acid motifs in a set of proteins.
#'
#' @param proteins A list of protein sequences read by `read.fasta` from the `seqinr` package.
#' @param Nmer An integer specifying the length of the amino acid motifs.
#' @param mode A character string, either "prop" for proportions or "counts" for raw counts.
#'
#' @return A data.table where each row represents a protein and columns reflect N-mer frequencies.
#' @export
#'
#' @examples
#' proteins <- read.fasta("~/Documents/TAIR10_pep_20101214_updated")
#' Nmer_protein_aa(proteins, Nmer = 2, mode = "counts")

Nmer_protein_aa <- function(proteins, Nmer, mode) {
  require(data.table)
  require(seqinr)

  Nmer_step <- Nmer - 1
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  motifs <- do.call(paste0, expand.grid(rep(list(amino_acids), Nmer)))

  total <- length(proteins)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  results <- rbindlist(lapply(1:total, function(i) {
    setTxtProgressBar(pb, i)
    seq <- toupper(paste0(unlist(proteins[[i]]), collapse = ""))

    walk <- sapply(Nmer:nchar(seq), function(j) {
      substr(seq, j - Nmer_step, j)
    })

    walk_table <- data.table(table(walk))

    if (mode == "prop") {
      walk_table$N <- prop.table(walk_table$N)
    } else if (mode == "counts") {
      walk_table$N <- walk_table$N
    } else {
      stop("'mode' must be either 'prop' or 'counts'")
    }

    walk_table[, row_id := 1]
    rotated_dt <- dcast(walk_table, row_id ~ walk, value.var = "N")
    rotated_dt[, row_id := NULL]
    rotated_dt[, Protein := names(proteins)[i]]

    return(rotated_dt)
  }), fill = TRUE)

  close(pb)

  results <- results[, c("Protein", setdiff(names(results), "Protein")), with=F]

  return(results)
}






