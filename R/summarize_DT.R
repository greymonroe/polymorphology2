#' Summarize Data Table
#'
#' This function summarizes the specified columns of a data.table by calculating
#' the mean, median, and standard error (SE) for each group, and also provides
#' the count of observations per group.
#'
#' @param DT A `data.table` object. The data to summarize.
#' @param SDcols A character vector of column names. These are the columns on
#'   which to compute summary statistics (mean, median, SE).
#' @param BYcols A character vector of column names. These columns will be used
#'   for grouping the data.
#'
#' @return A `data.table` containing the mean, median, SE, and count (`N`) for
#'   each group defined by `BYcols`.
#'
#' @details
#' The function computes the following statistics for each group:
#' - Mean (`_mean`)
#' - Median (`_median`)
#' - Standard Error (SE, `_se`), which is calculated as `sd / sqrt(n)`.
#' - Count of observations in each group (`N`).
#'
#' @examples
#' library(data.table)
#' DT <- data.table(id = c(1, 1, 2, 2, 3, 3),
#'                  value1 = c(5, 6, 7, 8, 9, 10),
#'                  value2 = c(50, 60, 70, 80, 90, 100))
#'
#' # Summarize by 'id' with 'value1' and 'value2' columns
#' summarize_DT(DT, SDcols = c("value1", "value2"), BYcols = "id")
#'
#' @export
summarize_DT <- function(DT, SDcols, BYcols) {
  DT[, c(setNames(lapply(.SD, mean), paste0(names(.SD), "_mean")),
         setNames(lapply(.SD, median), paste0(names(.SD), "_median")),
         setNames(lapply(.SD, function(x) sd(x) / sqrt(.N)), paste0(names(.SD), "_se")),
         .(N = .N)),  # Add the count per group
     by = BYcols, .SDcols = SDcols]
}
