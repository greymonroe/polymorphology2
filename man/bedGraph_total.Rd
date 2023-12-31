% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bedGraph_total.R
\name{bedGraph_total}
\alias{bedGraph_total}
\title{Calculate Total Depth of a BedGraph}
\usage{
bedGraph_total(bed)
}
\arguments{
\item{bed}{A data.table returned by the function read.bedGraph().
It must include the columns "CHROM", "START", "STOP", and "DEPTH".}
}
\value{
The total depth of the bedGraph.
}
\description{
This function calculates the total depth of a bedGraph.
The total depth is calculated by summing the product of DEPTH and the length of the interval (STOP - START) for each row.
}
\examples{
# Ensure bedGraph data.table is properly structured before using this function.
# Example:
# bed <- data.table(CHROM = c("chr1", "chr1"), START = c(1, 100), STOP = c(50, 150), DEPTH = c(2, 3))
# bedGraph_total(bed)

}
\seealso{
\code{\link[data.table]{data.table}}, \code{\link[data.table]{read.bedGraph}}
}
\keyword{genomics}
