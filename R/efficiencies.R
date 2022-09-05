#' @title Efficiencies
#' 
#' @description Extract the scores (optimal objective values) of the evaluated DMUs from a DEA / DEA fuzzy solution.
#' Note that these scores may not always be interpreted as efficiencies.
#' 
#' @param x DEA / DEA fuzzy object
#' @param ... ignored
#' 
#' @export 

efficiencies <- function(x,...) {
  UseMethod("efficiencies",x)
}

efficiencies.default <- function(x,...) {
  "Unknown class"
}
