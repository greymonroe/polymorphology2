#' N-mer Frequency
#'
#' Calculate the frequency of N-mers in given genomic features.
#'
#' @param features A `data.table` containing genomic features. It must include the columns "CHROM", "START", "STOP", and "ID".
#' @param fasta A named list representing the genome sequence. Names of the list should correspond to chromosome names, and each item in the list should be a character vector representing the sequence of the corresponding chromosome.
#' @param Nmer numeric value (2 would give dimers, 3 would give trimers, etc)
#' @param mode a string that can be either "prop" or "counts" to return either Nmer frequencies as total numbers or proportions.
#'
#' @return A `data.table` containing the frequency of each N-mer for each genomic feature.
#' @export
#'
#' @examples
#' Nmerfrequency(features, fasta, Nmer = 2, mode = "prop")

Nmerfrequency <- function(features, fasta, Nmer, mode = "prop") {

  required_cols <- c("CHROM", "START", "STOP", "ID")
  if (any(!(required_cols %in% colnames(features)))) {
    stop(paste("Features must contain the columns:", paste(required_cols, collapse = ", ")))
  }

  if(any(!(features$CHROM %in% names(fasta)))){
    warning(paste("The following CHROM from features are not found in names of fasta object:", paste(features$CHROM[!features$CHROM %in% names(fasta)], collapse = ", ")))
  }

  # Generate all possible N-mers
  motifs <- do.call(paste0, expand.grid(rep(list(c("A", "T", "C", "G")), Nmer)))

  # Create empty data.table
  empty_dt <- data.table(matrix(ncol = length(motifs), nrow = 0))
  setnames(empty_dt, motifs)

  Nmer_step=Nmer-1

  # Initialize a progress bar
  pb <- txtProgressBar(min = 0, max = nrow(features), style = 3)

  out <- rbindlist(lapply(seq_len(nrow(features)), function(i) {
    # Update the progress bar
    setTxtProgressBar(pb, i)

    x <- features[i]
    chrom <- x[["CHROM"]]
    start <- max(1, as.numeric(x[["START"]]))  # Ensuring start is at least 1
    stop  <- min(length(fasta[[chrom]]), as.numeric(x[["STOP"]]))  # Ensuring stop doesn't exceed the chromosome length
    seq <- toupper(paste0(unlist(fasta[[chrom]][start:stop]), collapse = ""))

    # allhits <- rbindlist(lapply(motifs, function(motif) {
    #   hits <- str_locate_all(seq, motif)
    #   data.table(motif, count = nrow(hits[[1]]))
    # }))
    #

    walk<-sapply(Nmer:nchar(seq), function(i){
      substr(seq, i-Nmer_step, i)
    })

    walk_table<-data.table(table(walk))
    rm("walk")
    if(mode == "prop"){
      walk_table$N <- prop.table(walk_table$N)
    } else if(mode == "counts"){
      walk_table$N <- walk_table$N
    } else {
      stop("'mode' must be either 'prop' or 'counts'")
    }

    walk_table[, row_id := 1]

    # Now we use dcast to "rotate" the data table
    rotated_dt <- dcast(walk_table, row_id ~ walk, value.var = "N")

    # Drop the row_id column
    rotated_dt[, row_id := NULL]

    rotated_dt[, ID := x[["ID"]]]

  }), fill = T)

  # Close the progress bar
  close(pb)

  out<-rbind(empty_dt, out, fill=T)
  out[ , names(out) := lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

  return(out)
}
