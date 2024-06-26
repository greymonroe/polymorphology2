% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tricontexts.R
\name{plot_tricontexts}
\alias{plot_tricontexts}
\title{Plots Trinucleotide Contexts of Mutations}
\usage{
plot_tricontexts(contexts, full = TRUE, trimer_freq = NULL)
}
\arguments{
\item{contexts}{A vector of mutation contexts. These contexts should be generated by the tricontexts() function.}

\item{full}{Logical flag indicating whether to plot the full context (TRUE) or only the mutations (FALSE).}

\item{trimer_freq}{data.table with columns TRI (phased trimer) and N (counts across region/genome)
Defaults to TRUE.}
}
\value{
A list containing two elements:
        context_table: A data.table object containing the mutation contexts and their respective counts.
        plot: A ggplot object of the mutation contexts.
}
\description{
This function plots the trinucleotide context for a set of mutation contexts.
The context is plotted either in a full mode or mutation only mode based on the full flag.
}
