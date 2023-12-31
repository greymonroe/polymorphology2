% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/features_chip_enrich.R
\name{features_chip_enrich}
\alias{features_chip_enrich}
\title{Calculate ChIP Enrichment for Genomic Features}
\usage{
features_chip_enrich(features, chipfile, inputfile)
}
\arguments{
\item{features}{A `data.table` containing genomic features. It must include the columns "CHROM", "START", "STOP", and "ID".}

\item{chipfile}{The name of the ChIP bedGraph file.}

\item{inputfile}{The name of the input bedGraph file.}
}
\value{
A `data.table` with the calculated enrichment for each feature.
}
\description{
This function calculates the ChIP enrichment for a set of genomic features using ChIP and input bedGraph files.
It reads the bedGraph files, calculates the overlaps with the features, and calculates the enrichment.
}
\examples{
# Ensure features data.table is properly structured before using this function.
# Example:
# features <- data.table(CHROM = c("chr1", "chr1", "chr2"), START = c(1, 100, 200), STOP = c(50, 150, 250), ID = c("feat1", "feat2", "feat3"))
# features_chip_enrich(features = features, chipfile = "chip.bedGraph", inputfile = "input.bedGraph")

}
\seealso{
\code{\link[data.table]{data.table}}, \code{\link[data.table]{read.bedGraph}}, \code{\link[data.table]{bedGraph_total}}, \code{\link[data.table]{features_in_features}}
}
\keyword{genomics}
