#' My Custom Theme
#'
#' This theme provides a transparent plot and panel background.
#'
#' @export
#'

my_theme <- function(base_size = 6, ...){
  theme_minimal(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.line = element_blank(),
      panel.border = element_blank(),
      ...
    )
}
