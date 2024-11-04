#' Box-Cox Transformation with Automatic Lambda Selection
#'
#' This function applies a Box-Cox transformation to a numeric vector, automatically selecting the optimal lambda
#' parameter to best approximate a normal distribution. The Box-Cox transformation is particularly useful for
#' stabilizing variance and making skewed data more normally distributed.
#'
#' @param x A numeric vector of data values to be transformed. All values must be positive. For data containing negative values, consider shifting the data to positive values before using this function.
#' @return A numeric vector of transformed values, using the optimal lambda for normality. The transformed values will approximate a normal distribution with stabilized variance.
#' @details The function first estimates the optimal lambda by maximizing the log-likelihood of the Box-Cox transformed values.
#' The optimal lambda is applied to transform the original data. If the optimal lambda is zero, a natural log transformation is used instead.
#'
#' @references
#' Box, G. E. P., & Cox, D. R. (1964). An Analysis of Transformations. \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 26(2), 211–243.
#'
#' @examples
#' # Example with a positively skewed dataset
#' set.seed(42)
#' skewed_data <- rexp(100, rate = 0.5)
#' # Apply the Box-Cox transformation
#' transformed_data <- boxcox_transform(skewed_data)
#' # View result
#' hist(transformed_data, main = "Box-Cox Transformed Data", col = "lightblue", border = "black")
#'
#' @importFrom MASS boxcox
#' @export
boxcox_transform <- function(x) {
  boxcox_result <- boxcox(lm(x ~ 1), lambda = seq(-5, 5, by = 0.1))
  optimal_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
  transformed_data <- if (optimal_lambda == 0) log(x) else (x ^ optimal_lambda - 1) / optimal_lambda
  return(transformed_data)
}
