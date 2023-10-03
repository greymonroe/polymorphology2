#' Format numbers with real scientific notation
#'
#' This function formats numeric values into a scientific notation with true exponents.
#' For example, instead of '1e-4', the output would render as \(1 \times 10^{-4}\) on a plot.
#'
#' @param x A numeric vector of values to be formatted.
#'
#' @return An expression vector with numbers formatted in scientific notation suitable for plot rendering.
#' @export
#'
#' @examples
#' real_sci_format(c(1000, 0.0001, 1e-7))
real_sci_format <- function(x) {
  format <- function(val) {
    if (is.na(val)) {
      return(as.character(NA))
    }
    if (val == 0) {
      return(expression(0))
    }
    e <- floor(log10(abs(val)))
    m <- val/10^e
    if (e == 0) {
      return(as.expression(m))
    }
    expr_text <- paste0(m, " %*% 10^", e)
    return(parse(text = expr_text))
  }
  sapply(x, format)
}
