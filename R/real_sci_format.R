#' Format numbers with real scientific notation
#'
#' This function formats numeric values into a scientific notation with true exponents.
#' For example, instead of '1e-4', the output would render as \(1 \times 10^{-4}\) on a plot.
#'
#' @param x A numeric vector of values to be formatted.
#' @param signif_N A numeric value for rounding small values to signif_N digits with signif, default is NULL
#'
#' @return An expression vector with numbers formatted in scientific notation suitable for plot rendering.
#' @export
#'
#' @examples
#' real_sci_format(c(1000, 0.0001, 1e-7))
real_sci_format <- function(x, signif_N=NULL) {
  format <- function(val) {
    if (is.na(val)) {
      return(as.character(NA))
    }
    if (val == 0) {
      return(expression(0))
    }
    if(!is.null(signif_N)) val<-signif(val, signif_N)
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
