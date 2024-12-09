#' Plot Grid with Adjustable Layout
#'
#' This function arranges a list of ggplot objects in a grid, either in columns or rows,
#' with adjustable relative sizes for each plot.
#'
#' @param plotlist A list of ggplot objects to be arranged.
#' @param type A character string specifying the layout type. Options are "cols" (side-by-side)
#'             or "rows" (stacked). Default is "cols".
#' @param relsize A numeric vector of the same length as \code{plotlist}, specifying the relative
#'                sizes of the plots. Default is equal sizing (\code{rep(1, length(plotlist))}).
#'
#' @return No explicit return value. The function directly draws the combined plot grid.
#'
#' @examples
#' library(ggplot2)
#'
#' # Example plots
#' p1 <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' p2 <- ggplot(mtcars, aes(gear, carb)) + geom_bar(stat = "identity")
#'
#' # Side-by-side with relative widths
#' plot_grid2(plotlist = list(p1, p2), type = "cols", relsize = c(5, 1))
#'
#' # Stacked with relative heights
#' plot_grid2(plotlist = list(p1, p2), type = "rows", relsize = c(2, 1))
#'
#' @export
plot_grid2 <- function(plotlist, type = "cols", relsize = rep(1, length(plotlist))) {
  if (!type %in% c("cols", "rows")) {
    stop("Invalid 'type'. Must be 'cols' or 'rows'.")
  }

  if (length(plotlist) != length(relsize)) {
    stop("Length of 'relsize' must match the number of plots in 'plotlist'.")
  }

  # Convert ggplot objects to grobs
  grobs <- lapply(plotlist, ggplotGrob)

  # Adjust widths or heights based on the type
  if (type == "cols") {
    for (i in seq_along(grobs)) {
      grobs[[i]]$widths <- grid::unit.pmax(grobs[[i]]$widths, grid::unit(relsize[i], "null"))
    }
    combined_grob <- do.call(cbind, grobs)
  } else if (type == "rows") {
    for (i in seq_along(grobs)) {
      grobs[[i]]$heights <- grid::unit.pmax(grobs[[i]]$heights, grid::unit(relsize[i], "null"))
    }
    combined_grob <- do.call(rbind, grobs)
  }

  # Draw the combined grob
  #grid::grid.newpage()
  grid::grid.draw(combined_grob)
}
