% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/INT.R
\name{INT}
\alias{INT}
\title{Rank-based Inverse Normal Transformation (INT)}
\usage{
INT(x)
}
\arguments{
\item{x}{A numeric vector of data values to be transformed.}
}
\value{
A numeric vector of transformed values, approximately normally distributed with mean 0 and standard deviation 1.
}
\description{
This function applies a rank-based inverse normal transformation (INT) to a vector of data, transforming it to follow an approximate normal distribution.
INT is often used to handle non-normal, skewed data by transforming it into a form suitable for statistical analyses requiring normality.
}
\details{
The function ranks the data, converts ranks to probabilities, and then applies the inverse normal cumulative distribution function.
INT is a common preprocessing step in genetics and other fields where traits may not follow a normal distribution.

This transformation technique is particularly useful in rare-variant association testing where phenotypic outliers and non-normality can affect the analysis (Auer, Reiner, & Leal, 2016).
}
\examples{
# Generate a skewed dataset
x <- rexp(100, rate = 0.5)
# Apply the INT transformation
transformed_x <- INT(x)
# View result
hist(transformed_x, main = "INT Transformed Data", col = "skyblue")
}
\references{
Auer, P., Reiner, A., & Leal, S. (2016). The effect of phenotypic outliers and non-normality on rare-variant association testing. \emph{Eur J Hum Genet, 24}, 1188–1194. https://doi.org/10.1038/ejhg.2015.270
}
