#' @title deadata_fuzzy class check.
#'   
#' @description Checks whether an R object is of deadata_fuzzy class or not.
#' @usage is.deadata_fuzzy(x)
#' 
#' @param x Any \bold{R} object.
#' @return Returns \code{TRUE} if its argument is a deadata_fuzzy object (that is, has "deadata_fuzzy" 
#'   amongst its classes) and \code{FALSE} otherwise.
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
#' @export


is.deadata_fuzzy <- function(x) {
  
  inherits(x, "deadata_fuzzy")
  
}