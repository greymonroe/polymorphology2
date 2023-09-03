#' plot_bins Function
#'
#' This function takes a data.table, two column names, a number of bins, and an x-axis type as input.
#' It then calculates the mean and standard error of the y-variable for each bin of the x-variable.
#' Finally, it returns a ggplot object showing the binned data.
#'
#' @param data A data.table containing the data.
#' @param yvar Character name of the column in data to be used as the y-variable.
#' @param xvar Character name of the column in data to be used as the x-variable.
#' @param bins Number of bins to divide the x-variable into.
#' @param xaxis Type of x-axis to be used in the plot. Either "bins" or "numeric".
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
#' result <- plot_bins(dt, "g1001", "enrich", 10, "bins")
#' print(result$plot)
#' }
#'
plot_bins <- function(data, yvar, xvar, bins, xaxis = "numeric") {
  require(ggplot2)
  require(data.table)

  bins_data <- data[is.finite(get(yvar)), .(y = mean(get(yvar), na.rm=T), y_se = sd(get(yvar), na.rm=T) / sqrt(.N), N = .N),
                    by = .(x_cut = cut2(get(xvar), g = bins))]
  bins_data$x_num <- as.numeric(bins_data$x_cut)

  if (xaxis == "numeric") {
    plot <- ggplot(bins_data, aes(x = x_num, y = y)) +
      geom_point() +
      geom_errorbar(aes(ymin = y - 2 * y_se, ymax = y + 2 * y_se), width = 0) +
      scale_x_continuous(name = xvar) +
      scale_y_continuous(name = yvar)+
      theme_classic()
  } else if (xaxis == "bins") {
    plot <- ggplot(bins_data, aes(x = x_cut, y = y)) +
      geom_point() +
      theme_classic()+
      geom_errorbar(aes(ymin = y - 2 * y_se, ymax = y + 2 * y_se), width = 0) +
      scale_y_continuous(name = yvar) +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  } else {
    stop("Invalid value for xaxis. Choose either 'numeric' or 'bins'.")
  }

  return(list(bins_data = bins_data, plot = plot))
}
