#' plot_bins Function
#'
#' This function takes a data.table, two column names, and a number of bins as input.
#' It then calculates the mean and standard error of the y-variable for each bin of the x-variable.
#' Finally, it returns a ggplot object showing the binned data.
#'
#' @param data A data.table containing the data.
#' @param yvar Character name of the column in data to be used as the y-variable.
#' @param xvar Character name of the column in data to be used as the x-variable.
#' @param bins Number of bins to divide the x-variable into.
#'
#' @return A list containing the binned data and the ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar scale_x_continuous scale_y_continuous
#' @importFrom data.table .N
#' @importFrom Hmisc cut2
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(enrich = rnorm(100), g1001 = rnorm(100))
#' result <- plot_bins(dt, "g1001", "enrich", 10)
#' print(result[[2]])
#' }
#'
plot_bins <- function(data, yvar, xvar, bins) {
  require(ggplot2)
  require(data.table)

  bins_data <- data[, .(y = mean(get(yvar)), y_se = sd(get(yvar)) / sqrt(.N), N = .N),
                    by = .(x_cut = cut2(get(xvar), g = bins))]
  bins_data$x_num <- as.numeric(bins_data$x_cut)

  plot <- ggplot(bins_data, aes(x = x_num, y = y)) +
    geom_point() +
    geom_errorbar(aes(ymin = y - 2 * y_se, ymax = y + 2 * y_se), width = 0) +
    scale_x_continuous(name = xvar) +
    scale_y_continuous(name = yvar)

  return(list(bins_data=bins_data, plot=plot))
}
