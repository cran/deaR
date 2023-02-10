#' @title deadata class print method
#'   
#' @description Print method for \code{deadata} class.
#' 
#' @param x A \code{deadata} object (as returned by \code{make_deadata} function).
#' @param ... For compatibility issues.
#'   
#' @author 
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#' 
#' \strong{Vicente Bolós} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benítez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#' 
#' @method print deadata
#' @export


print.deadata <- function(x, ...){
  print.default(x)
}
