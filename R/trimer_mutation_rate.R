#' Trimer Mutation Rate
#'
#' Calculate the mutation rate for each trimer.
#'
#' @param tricontexts A vector generated from the `tricontexts()` function.
#' @param trimer_counts A `data.table` generated from the `trimerfrequencycomp(Nmerfrequency(Nmer=3, mode="counts"))` function.
#'
#' @return A `data.table` containing the mutation rate for each trimer.
#' @export
#'
#' @examples
#' trimer_mutation_rate(tricontexts, trimer_counts)

trimer_mutation_rate <- function(tricontexts, trimer_counts) {

  # Create a table of the count of each trimer in the 'tricontexts' vector
  trimut_count <- data.table((table(trimer = substr(tricontexts, 1, 3))))
  colnames(trimut_count)[2] <- "mutation_counts"

  # Calculate the proportion of each trimer
  trimut_count$mutation_prop <- prop.table(trimut_count$mutation_counts)

  # Create a table of the count of each trimer in the 'trimer_counts' data.table
  tri_counts <- data.table(trimer = colnames(trimer_counts)[2:33], trimer_count = colSums(trimer_counts[, 2:33]))

  # Calculate the proportion of each trimer
  tri_counts$trimer_prop <- prop.table(tri_counts$trimer_count)

  # Merge the two data.tables by 'trimer'
  tri_counts <- merge(tri_counts, trimut_count)

  # Calculate the expected count of each trimer
  tri_counts$expected <- chisq.test(tri_counts$mutation_counts, p = prop.table(tri_counts$trimer_count))$expected

  # Calculate the observed/expected ratio
  tri_counts$oe <- tri_counts$mutation_counts / tri_counts$expected

  # Calculate the mutation rate
  tri_counts$rate <- tri_counts$mutation_counts / tri_counts$trimer_count

  return(tri_counts)
}
