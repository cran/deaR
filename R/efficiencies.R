#' @title Objective scores
#' 
#' @description Extract the scores (optimal objective values) of the evaluated
#' DMUs from a conventional, fuzzy or stochastic DEA solution. Note that these
#' scores may not always be interpreted as efficiencies.
#' 
#' @param x An object of class \code{dea}, \code{dea_fuzzy} or \code{dea_stoch}.
#' @param ... ignored.
#' 
#' @export 

efficiencies <- function(x, ...) {
  UseMethod("efficiencies", x)
}

#' @export
#' 
efficiencies.default <- function(x, ...) {
  "Unknown class"
}
