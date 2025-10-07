#' @title read_malmquist
#'  
#' @description This function is deprecated. Use \code{make_malmquist} instead.
#' 
#' @usage read_malmquist(datadea,
#'                nper = NULL,
#'                percol = NULL,
#'                arrangement  = c("horizontal", "vertical"),
#'                ...) 
#'              
#' @param datadea Data frame with DEA data.
#' @param nper Number of time periods (with dataset in wide format).
#' @param percol Column of time period (with dataset in long format).
#' @param arrangement Horizontal with data in wide format. Vertical with data in long format.
#' @param ... Other options to be passed to the \code{make_deadata} function.
#' 
#' @returns An object of class \code{deadata}. Function \code{malmquist_index} can be applied.
#'
#' @export

read_malmquist <- function(datadea,
                           nper = NULL,
                           percol = NULL,
                           arrangement  = c("horizontal","vertical"),
                           ...) {
  
  .Deprecated("make_malmquist")
  
  make_malmquist(datadea = datadea,
                 nper = nper,
                 percol = percol,
                 arrangement  = arrangement,
                 ...)

}