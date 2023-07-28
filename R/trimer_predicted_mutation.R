#' Predicted Mutation Rate for Trimers
#'
#' This function calculates the expected mutation count and total for each trimer.
#'
#' @param trimer_counts A `data.table` generated from the `trimerfrequencycomp(Nmerfrequency(Nmer=3, mode="counts"))` function.
#' @param trimer_mutation_rates A `data.table` from the `trimer_mutation_rate()` function.
#'
#' @return A `data.table` containing the expected mutation count and total for each trimer.
#' @export
#'
#' @examples
#' trimer_predicted_mutation(trimer_counts, trimer_mutation_rates)
trimer_predicted_mutation <- function(trimer_counts, trimer_mutation_rates) {

  # Ensure trimer names from trimer_counts and trimer_mutation_rates are the same
  if(!all(colnames(trimer_counts)[2:ncol(trimer_counts)] %in% trimer_mutation_rates$trimer)) {
    stop("Trimers in trimer_counts do not match with trimers in trimer_mutation_rates")
  }

  # Calculate the expected mutation counts for each trimer
  trimer_counts$counts_exp <- apply(trimer_counts[, trimer_mutation_rates$trimer, with = F], 1, function(x){
    expected <- sum(x * trimer_mutation_rates$rate)
    if(is.na(expected)) {
      warning(paste("Expected count for trimer", names(x), "could not be calculated. Check if 'rate' is defined for this trimer in 'trimer_mutation_rates'"))
    }
    return(expected)
  })

  # Calculate the total for each trimer
  trimer_counts$trimer_total <- apply(trimer_counts[, trimer_mutation_rates$trimer, with = F], 1, function(x){
    return(sum(x))
  })

  return(trimer_counts)
}
