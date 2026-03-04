#' QQ Plot for P-values
#'
#' Generates a QQ plot to visualize the distribution of observed versus expected -log10 transformed p-values.
#'
#' @param p_values A numeric vector of p-values to be plotted.
#'
#' @return A ggplot object representing the QQ plot.
#'
#' @details This function calculates the expected and observed -log10(p-values) and
#' plots them on a QQ plot. The diagonal line represents the null hypothesis of
#' no association (y = x). Points deviating from the line suggest potential associations.
#'
#' @examples
#' # Example usage:
#' p_values <- runif(100, min = 1e-8, max = 1)
#' qqplot(p_values)
#'
#' @export
#'
qqplot<-function(p_values){

  # Calculate expected and observed -log10(p-values)
  oe<-data.table(p_values=sort(p_values))
  oe$observed<- -log10((oe$p_values))
  oe$expected <- -log10(ppoints(length(oe$observed)))


  # Create the QQ plot
  qqp<-ggplot(data =oe , aes(x = expected, y = observed)) +
    geom_point(size = 0.01, col="red") +
    scale_color_manual(values=c("gray","red"), guide="none")+
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.25) + # Diagonal line (y = x)
    labs(x = expression(Expected ~ -log[10](p)), y = expression(Observed ~ -log[10](p))) +
    theme_void(base_size = 6)+
    coord_fixed()+
    theme(plot.title = element_text(size=6), axis.line = element_line(), axis.title = element_blank(), panel.grid = element_blank())
  return(qqp)

}

