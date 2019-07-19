#' @title Efficiencies
#' @description Extracts efficiencies from  dea/dea_fuzzy objects.
#' @param x dea / deafuzzy object
#' @param ... ignored
#' @export 
efficiencies <- function(x,...) {
  UseMethod("efficiencies",x)
}

efficiencies.default <- function(x,...) {
  "Unknown class"
}
