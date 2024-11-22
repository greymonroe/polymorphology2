#' N-mer Frequency
#'
#' Calculate the frequency of N-mers in given genomic features.
#'
#' @param features A `data.table` containing genomic features. It must include the columns "CHROM", "START", "STOP", and "ID".
#' @param fasta A named list representing the genome sequence. Names of the list should correspond to chromosome names, and each item in the list should be a character vector representing the sequence of the corresponding chromosome.
#' @param Nmer numeric value (2 would give dimers, 3 would give trimers, etc)
#' @param mode a string that can be either "prop" or "counts" to return either Nmer frequencies as total numbers or proportions.
#' @param chunk_size an integer, the number of rows per chunk for faster processing and lower memory use
#'
#' @return A `data.table` containing the frequency of each N-mer for each genomic feature.
#' @export


Nmerfrequency <- function(features, fasta, Nmer, mode = "prop", parallel = TRUE, chunk_size = 500) {
  required_cols <- c("CHROM", "START", "STOP", "ID")
  if (any(!(required_cols %in% colnames(features)))) {
    stop(paste("Features must contain the columns:", paste(required_cols, collapse = ", ")))
  }
  if (any(!(features$CHROM %in% names(fasta)))) {
    warning(paste("The following CHROM from features are not found in names of fasta object:",
                  paste(features$CHROM[!features$CHROM %in% names(fasta)], collapse = ", ")))
  }

  # Generate all possible N-mers
  motifs <- do.call(paste0, expand.grid(rep(list(c("A", "T", "C", "G")), Nmer)))
  empty_dt <- data.table(matrix(0, ncol = length(motifs), nrow = 0))
  setnames(empty_dt, motifs)

  Nmer_step <- Nmer - 1

  getNmercounts <- function(i, subfeatures) {
    x <- subfeatures[i]
    chrom <- x[["CHROM"]]
    start <- as.numeric(x[["START"]])
    stop <- as.numeric(x[["STOP"]])

    seq <- toupper(fasta[[chrom]][start:stop])  # Faster substring extraction
    if (Nmer > 1) {
      seq <- toupper(paste0(unlist(fasta[[chrom]][start:stop]), collapse = ""))
      walk<-sapply(Nmer:nchar(seq), function(i){
        substr(seq, i-Nmer_step, i)
      }) } else {
        walk <- seq
      }

    walk_table <- as.data.table(table(walk))
    if (mode == "prop") {
      walk_table$N <- prop.table(walk_table$N)
    }

    walk_table[, row_id := 1]
    rotated_dt <- dcast(walk_table, row_id ~ walk, value.var = "N")
    rotated_dt[, row_id := NULL]
    rotated_dt[, ID := x[["ID"]]]

    return(rotated_dt)
  }

  # Determine the number of chunks
  num_chunks <- ceiling(nrow(features) / chunk_size)
  out <- rbindlist(pblapply(seq_len(num_chunks), function(chunk) {
    # Define the range of rows for this chunk
    start_row <- (chunk - 1) * chunk_size + 1
    end_row <- min(chunk * chunk_size, nrow(features))
    numrows<-length(start_row:end_row)
    subfeatures<-features[start_row:end_row]
    # Process this chunk in parallel if parallel = TRUE
    chunk_results <- if (parallel) {
      parallel::mclapply(1:numrows, function(i) {
        getNmercounts(i, subfeatures)
      }, mc.cores = detectCores() -1)
    } else {
      lapply(1:numrows, function(i) {
        getNmercounts(i, subfeatures)
      })
    }

    # Combine results from this chunk
    rbindlist(chunk_results, fill = TRUE)
  }), fill = TRUE)

  # Add empty columns to ensure all N-mers are represented
  out <- rbind(empty_dt, out, fill = TRUE)
  out[, (names(out)) := lapply(.SD, function(x) ifelse(is.na(x), 0, x))]

  return(out)
}
