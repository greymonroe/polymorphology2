#' Plot Variable in Relation to Features
#'
#' This function generates a plot of a specified variable in relation to features in a set of genomic feature windows.
#' It allows for calculation and plotting of the mean or percentage of a variable.
#'
#' @param feature_windows A `data.table` object created by the `feature_windows` function. It should include the columns "REGION", "RELATIVEPOS", and "LENGTH".
#' @param variable A string specifying the name of the column in the `feature_windows` data.table to be plotted.
#' @param mode A string specifying the mode of calculation: "mean" for mean of the variable, "percent" for the sum of the variable values divided by the total length.
#'
#' @return A `ggplot` object representing the plot of the variable in relation to features.
#'
#' @examples
#' # Ensure 'feature_windows' data.table is properly structured before using this function.
#' # Example:
#' # feature_windows <- feature_windows(features = features, chipfile = chipfile, inputfile = inputfile, mode = "mean", value="DEPTH")
#' # plot_feature_windows(feature_windows = feature_windows, variable = "DEPTH", mode = "mean")
#'
#' @export
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_line}}, \code{\link[ggplot2]{geom_vline}}, \code{\link[ggplot2]{theme_classic}}, \code{\link[ggplot2]{scale_x_continuous}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline theme_classic scale_x_continuous
#'
#' @keywords genomics, plotting

plot_feature_windows <- function(feature_windows, variable, mode) {

  if(!(mode %in% c("mean", "percent"))) {
    stop("Invalid 'mode'. Must be one of 'mean', 'percent'.")
  }

  if(!(variable %in% colnames(feature_windows))) {
    stop("The specified 'variable' column does not exist in the 'feature_windows' data.table.")
  }

  if (mode == "mean") {
    summary <- feature_windows[, .(y = mean(get(variable), na.rm=T)), by = .(REGION, RELATIVEPOS)]
  } else if (mode == "percent") {
    summary <- feature_windows[, .(y = sum(get(variable)) / sum(LENGTH)), by = .(REGION, RELATIVEPOS)]
  }

  maxpos <- max(summary$RELATIVEPOS)

  plot <- ggplot(summary, aes(x = RELATIVEPOS, y = y)) +
    geom_line(col = "red") +
    geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2), linetype = "dashed") +
    theme_classic(base_size = 6) +
    scale_x_continuous(breaks = c(1, maxpos / 3, maxpos / 3 * 2, maxpos), labels = c("-2kb", "START", "STOP", "+2kb"))

  return(list(summary=summary, plot=plot))
}
