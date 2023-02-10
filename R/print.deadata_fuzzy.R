#' @title deadata_fuzzy class print method
#'   
#' @description Print method for \code{deadata_fuzzy} class.
#' 
#' @param x A \code{deadata_fuzzy} object (as returned by \code{make_deadata_fuzzy} function).
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
#' @method print deadata_fuzzy
#' @export


print.deadata_fuzzy <- function(x, ...){
  print.default(x)
}
