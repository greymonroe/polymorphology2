% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trimer_predicted_mutation.R
\name{trimer_predicted_mutation}
\alias{trimer_predicted_mutation}
\title{Predicted Mutation Rate for Trimers}
\usage{
trimer_predicted_mutation(trimer_counts, trimer_mutation_rates)
}
\arguments{
\item{trimer_counts}{A `data.table` generated from the `trimerfrequencycomp(Nmerfrequency(Nmer=3, mode="counts"))` function.}

\item{trimer_mutation_rates}{A `data.table` from the `trimer_mutation_rate()` function.}
}
\value{
A `data.table` containing the expected mutation count and total for each trimer.
}
\description{
This function calculates the expected mutation count and total for each trimer.
}
\examples{
trimer_predicted_mutation(trimer_counts, trimer_mutation_rates)
}
