% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonsyn_syn_exp3.R
\name{nonsyn_syn_exp3}
\alias{nonsyn_syn_exp3}
\title{Nonsynonymous to Synonymous Mutation Expectation}
\usage{
nonsyn_syn_exp3(muts, genome, CDS, reps)
}
\arguments{
\item{muts}{A data.table with REF and ALT columns indicating the reference and alternate nucleotides.}

\item{genome}{A fasta object read in by seqinr read.fasta()}

\item{CDS}{A coding sequence fasta object read by seqinr.}

\item{reps}{Number of repetitions for the simulation.}
}
\value{
A list containing the results for each repetition and the overall ratio of nonsynonymous to synonymous mutations.
}
\description{
This function calculates the expected ratio of nonsynonymous to synonymous mutations
based on a given mutation dataset and a coding sequence, accounting for base pair composition of genome. For pre-calculated bp_freq for a genome sequence use nonsyn_syn_exp4()
}
