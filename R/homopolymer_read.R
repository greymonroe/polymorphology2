#' Identify Homopolymers in a Sequence
#'
#' This function identifies runs of homopolymers in a given sequence. It uses
#' run length encoding (RLE) to identify and count consecutive occurrences
#' of the same base pair (BP) in the sequence.
#'
#' @param seq A vector representing a sequence in which to identify homopolymers.
#' @param size The minimum length of homopolymers to identify.
#'
#' @return A data.table of identified homopolymers that have a length greater than
#'         or equal to the specified size. Each row represents a homopolymer,
#'         with columns for the base pair (var), start position (START), stop
#'         position (STOP), and length (LENGTH).
#'
#' @examples
#' # Define a vector representing a sequence
#' seq <- c("A", "A", "A", "T", "G", "G", "G", "G")
#' # Call the function to identify homopolymers
#' homopolymers <- homopolymer_read(seq, 2)
#'
#' @export
#'
#' @seealso \code{\link[rle]{rle}}
#'
#' @keywords genomics


homopolymer_read <- function(seq, size){
  # Calculate run length encoding
  s <- rle(seq)

  # Calculate stop positions
  v <- cumsum(s$lengths)

  # Construct data.table of runs
  runs <- data.table('BP' = toupper(s$values), 'START' = v + 1 - s$lengths, 'STOP' = v)

  # Calculate length of each run
  runs$LENGTH <- runs$STOP - runs$START + 1

  # Filter runs to get homopolymers of the specified size or larger
  homopolymers <- runs[LENGTH >= size]

  return(homopolymers)
}

