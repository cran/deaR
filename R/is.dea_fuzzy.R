#' @title dea_fuzzy class check.
#'   
#' @description Checks whether an R object is of dea_fuzzy class or not.
#' @usage is.dea_fuzzy(x)
#' 
#' @param x Any \bold{R} object.
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
#' @returns Returns \code{TRUE} if its argument is a \code{dea_fuzzy} object (that is, has "dea_fuzzy" 
#'   amongst its classes) and \code{FALSE} otherwise.
#'   
#' @export


is.dea_fuzzy <- function(x) {
  
  inherits(x, "dea_fuzzy")
  
}