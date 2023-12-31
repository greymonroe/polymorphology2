% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_windows.R
\name{feature_windows}
\alias{feature_windows}
\title{Create windows around genome features}
\usage{
feature_windows(features, breaks, dist, directed, IDcol)
}
\arguments{
\item{features}{A data.table with CHROM START and STOP columns. If directed=T,
it must contain a DIRECTION column.}

\item{breaks}{The number of windows to break the feature into.}

\item{dist}{distance from start and stop to make windows upstream and downstream}

\item{directed}{Logical. If TRUE, directionality of features is considered.}

\item{IDcol}{Column name in features data.table used to identify the feature.}
}
\value{
A data.table of windows with their respective features.
}
\description{
This function creates windows around the genome features. The regions upstream
and downstream of the features are split into windows of size 2000/breaks. If
directed = T, the output windows consider the direction of the features.
}
\examples{
# Ensure features is a data.table with proper structure before using this function.
}
