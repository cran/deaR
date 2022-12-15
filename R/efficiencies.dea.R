#' @title Objective scores
#'   
#' @description Extract the scores (optimal objective values) of the evaluated
#' DMUs from a conventional DEA solution. Note that these scores may not always be
#' interpreted as efficiencies.
#' 
#' @param x Object of class \code{dea} obtained with some of the conventional DEA
#' \code{model_*} functions.
#' @param ... Other options (for compatibility reasons).
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
#' @references 
#' Tomkins, C.; Green, R. (1988). “An Experiment in the Use of Data Envelopment
#' Analysis for Evaluating the Efficiency of UK University Departments of
#' Accounting”. Financial Accountability and Management 4(2): 147.
#' 
#' @examples 
#' # Replication results model DEA1 in Tomkins and Green (1988)
#' data("Departments")
#' # Calculate Total income
#' Departments$Total_income <- Departments[, 5] + Departments[, 6] + Departments[, 7] 
#' data_DEA1 <- make_deadata(Departments,
#'                           inputs = 9,
#'                           outputs = c(2, 3, 4, 12))
#' result <- model_basic(data_DEA1,
#'                       orientation = "io",
#'                       rts = "crs")
#' efficiencies(result) # Table 3 (p.156) 
#'  
#' @export

efficiencies.dea <-
  function(x, ...) {
    deasol <- x
    if ("efficiency" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$efficiency))
      if (length(deasol$DMU[[1]]$efficiency) > 1) {
        eff  <- do.call(rbind, lapply(deasol$DMU, function(x)
          x$efficiency))
        mean_eff <- unlist(lapply(deasol$DMU, function(x)
          x$mean_efficiency))
        eff <- cbind(eff, mean_eff)
      }
    } else if ("beta" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$beta))
    } else if ("delta" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$delta))
    } else if ("gamma" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$gamma))
    } else if ("objval" %in% names(deasol$DMU[[1]])) {
      eff <- unlist(lapply(deasol$DMU, function(x)
        x$objval))
    } else {
      stop("No efficiency/beta/delta/objval parameters in this solution!")
    }
    
    return(round(eff, 6))
    
}