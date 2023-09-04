#' Plot Random Forest Feature Importance
#'
#' This function takes a random forest model object and plots the feature importance
#' as a horizontal bar plot, sorted from most to least important.
#'
#' @param rf_model A random forest model object, typically created using the `randomForest` function.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{importance_data}: A data.table containing the importance scores for each feature.
#'   \item \code{plot}: A ggplot2 object representing the feature importance plot.
#' }
#'
#' @examples
#' \dontrun{
#' library(randomForest)
#' data(iris)
#' rf <- randomForest(Species ~ ., data = iris, importance = TRUE)
#' result <- plot_rf_importance(rf)
#' print(result$plot)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal
#' @importFrom randomForest importance
#' @importFrom data.table data.table
plot_rf_importance <- function(rf_model) {
  # Extract importance scores
  importance_data <- as.data.frame(importance(rf_model))

  # Add a column for variable names
  importance_data$Variable <- rownames(importance_data)

  # Sort by importance
  importance_data <- importance_data[order(importance_data$MeanDecreaseGini, decreasing = TRUE), ]

  # Create a horizontal bar plot
  p <- ggplot(importance_data, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Feature", y = "Importance Score", title = "Random Forest Feature Importance") +
    theme_minimal()

  return(list(importance_data=data.table(importance_data), plot=p))
}
