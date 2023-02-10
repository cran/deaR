#' @title Lambdas
#'   
#' @description Extract the lambdas of the DMUs from a dea or dea_fuzzy solution.
#' 
#' @usage lambdas(deasol)
#' 
#' @param deasol Object of class \code{dea} or \code{dea_fuzzy} obtained with
#' some of the DEA model functions.
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
#' @examples 
#' data("Coll_Blasco_2006")
#' data_example <- make_deadata(Coll_Blasco_2006,
#'                              ni = 2, 
#'                              no = 2)
#' result <- model_multiplier(data_example, 
#'                            orientation = "io",
#'                            rts = "crs")
#' lambdas(result)
#'  
#' @export

lambdas <-
  function(deasol) {
  
  if (is.dea(deasol)) {
    
    if ("lambda" %in% names(deasol$DMU[[1]])) {
      lamb  <- do.call(rbind, lapply(deasol$DMU, function(x)
        x$lambda))
      return(round(lamb, 5))
    } else {
      stop("No lambda parameters in this solution!")
    }

  } else if (is.dea_fuzzy(deasol)) {
    
    dmunames_eval <- names(deasol$dmu_eval)
    dmunames_ref <- names(deasol$dmu_ref)
    nde <- length(deasol$dmu_eval)
    ndr <- length(deasol$dmu_ref)
    
    if (grepl("kaoliu", deasol$modelname)) {
      
      nalpha <- length(deasol$alpha)
      
      if ("lambda" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) {
        
        lamb.W <- array(0,
                        dim = c(nde, ndr, nalpha),
                        dimnames = list(dmunames_eval, dmunames_ref, names(deasol$alphacut)))
        lamb.B <- lamb.W
        
        for (i in 1:nalpha) {
          lamb.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
            x$lambda))
          lamb.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
            x$lambda))
        }
        return(list(Worst = lamb.W, Best = lamb.B))
        
      } else {
        stop("No lambda parameters in this solution!")
      }
      
    } else if (grepl("possibilistic", deasol$modelname)) {
      
      nh <- length(deasol$h)
      
      if ("lambda" %in% names(deasol$hlevel[[1]]$DMU[[1]])) {
        
        lamb <- array(0,
                      dim = c(nde, ndr, nh),
                      dimnames = list(dmunames_eval, dmunames_ref, names(deasol$hlevel)))
        
        for (i in 1:nh) {
          lamb[, , i] <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x)
            x$lambda))
        }
        return(round(lamb, 5))
        
      } else {
        stop("No lambda parameters in this solution!")
      }
      
    }else{
      stop("Guo-Tanaka model does not have lambdas!")
      
    }
    
  } else {
    
    stop("Input should be a dea or dea_fuzzy class object!")
    
  }
  
}