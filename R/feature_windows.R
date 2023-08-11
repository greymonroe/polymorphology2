#' Create windows around genome features
#'
#' This function creates windows around the genome features. The regions upstream
#' and downstream of the features are split into windows of size 2000/breaks. If
#' directed = T, the output windows consider the direction of the features.
#'
#' @param features A data.table with CHROM START and STOP columns. If directed=T,
#' it must contain a DIRECTION column.
#' @param breaks The number of windows to break the feature into.
#' @param directed Logical. If TRUE, directionality of features is considered.
#' @param IDcol Column name in features data.table used to identify the feature.
#'
#' @return A data.table of windows with their respective features.
#'
#' @export
#' @examples
#' # Ensure features is a data.table with proper structure before using this function.
feature_windows <- function(features, breaks, dist, directed, IDcol){

  # Check if features is a data.table
  if (!is.data.table(features)){
    stop("Error: 'features' should be a data.table")
  }

  # Check if necessary columns exist in features
  required_columns <- c("CHROM", "START", "STOP", IDcol)
  if(directed){
    required_columns <- c(required_columns, "DIRECTION")
  }
  if(any(!(required_columns %in% names(features)))){
    stop(paste("features object needs to have", toString(required_columns), "columns"))
  }

  # Check if START and STOP columns are numeric
  if(!is.numeric(features[["START"]]) | !is.numeric(features[["STOP"]])){
    stop("Error: START and STOP columns should be numeric")
  }

  # Check if DIRECTION column only contains "+" and "-" when directed = T
  if(directed && !all(features[["DIRECTION"]] %in% c("+", "-"))){
    stop("Error: When directed = T, DIRECTION column should only contain '+' and '-' values")
  }

  # Create progress bar
  #pb <- txtProgressBar(min = 0, max = nrow(features), style = 3)

  windows <- rbindlist(apply(features, 1, function(x) {

    # Update progress bar
    #setTxtProgressBar(pb, which(features[[IDcol]] == x[IDcol]))

    chrom = x["CHROM"]
    body_starts = round(seq(as.numeric(x["START"]), as.numeric(x["STOP"]), length.out=breaks+1)[-(breaks+1)])
    body_stops = round(seq(as.numeric(x["START"]), as.numeric(x["STOP"]), length.out=breaks+1)[-1])
    upstream_starts = seq(as.numeric(x["START"])-dist, as.numeric(x["START"]), length.out=breaks+1)[-(breaks+1)]
    upstream_stops = seq(as.numeric(x["START"])-dist, as.numeric(x["START"]), length.out=breaks+1)[-1]
    downstream_starts = seq(as.numeric(x["STOP"]), as.numeric(x["STOP"])+dist, length.out=breaks+1)[-(breaks+1)]
    downstream_stops = seq(as.numeric(x["STOP"]), as.numeric(x["STOP"])+dist, length.out=breaks+1)[-1]

    out = data.table(
      CHROM = x["CHROM"],
      START = c(upstream_starts, body_starts, downstream_starts),
      STOP = c(upstream_stops, body_stops, downstream_stops),
      REGION = c(rep("upstream", length(upstream_starts)),
                 rep("gene body", length(body_starts)),
                 rep("downstream", length(downstream_starts)))
    )

    out[, (IDcol) := x[IDcol]]

    out$RELATIVEPOS <- 1:nrow(out)
    out$LENGTH <- out$STOP - out$START

    if(directed == T){
      direction = x["DIRECTION"]
      if(direction=="-"){
        out$RELATIVEPOS <- rev(out$RELATIVEPOS)
        out$REGION <- rev(out$REGION)
      }
    }

    return(out)

  }), fill=TRUE)

  # Close progress bar
  #close(pb)

  # Add 'ID' column
  windows[, ID := 1:.N]

  setkey(windows, CHROM, START, STOP)
  return(windows)
}
