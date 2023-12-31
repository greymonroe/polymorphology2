% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homopolymer_read.R
\name{homopolymer_read}
\alias{homopolymer_read}
\title{Identify Homopolymers in a Sequence}
\usage{
homopolymer_read(seq, size)
}
\arguments{
\item{seq}{A vector representing a sequence in which to identify homopolymers.}

\item{size}{The minimum length of homopolymers to identify.}
}
\value{
A data.table of identified homopolymers that have a length greater than
        or equal to the specified size. Each row represents a homopolymer,
        with columns for the base pair (var), start position (START), stop
        position (STOP), and length (LENGTH).
}
\description{
This function identifies runs of homopolymers in a given sequence. It uses
run length encoding (RLE) to identify and count consecutive occurrences
of the same base pair (BP) in the sequence.
}
\examples{
# Define a vector representing a sequence
seq <- c("A", "A", "A", "T", "G", "G", "G", "G")
# Call the function to identify homopolymers
homopolymers <- homopolymer_read(seq, 2)

}
\seealso{
\code{\link[rle]{rle}}
}
\keyword{genomics}
