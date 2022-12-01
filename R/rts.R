#' @title Returns to scale
#'   
#' @description Extract the returns to scale. 
#'
#' @usage rts(deamodel,
#'            thr = 1e-4)
#'            
#' @param deamodel Object of class dea obtained with some of the dea functions.
#' @param thr Threshold for the tolerance for considering something = 1.
#' Defaults to 1e-4.
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
#' result <- model_basic(data_example, 
#'                       orientation = "io", 
#'                       rts ="crs")
#' rts(result)
#'  
#' @export

rts <- function(deamodel,
                thr =  1e-4) {
  
  if (!is.dea(deamodel)) {
    stop("Input should be a dea class object!")
  }
  rts <- NULL
  if (!deamodel$modelname %in% c("multiplier")) {
    lamb <- lambdas(deamodel)
    lambsum <- rowSums(lamb)
    if (deamodel$orientation %in% c("io","oo")) {
      switch(deamodel$rts,
             crs = {
               rts <- ifelse(lambsum > 1 + thr , "Decreasing",
                             ifelse(abs(lambsum - 1) < thr, "Constant", "Increasing"))
             },
             vrs = {
               rts <- ifelse(abs(lambsum - 1) < thr, "Variable","Variable")
             },
             nirs = {
               rts <- ifelse(lambsum < 1 - thr, "Decreasing","Constant")
             },
             ndrs = {
               rts <- ifelse(lambsum > 1 + thr, "Increasing","Constant")
             },
             warning("RTS with General returns to scale are not implemented yet!")
             
      )
      if (deamodel$rts != "grs") {
        res <- data.frame(lambsum = lambsum, rts = rts)
      } else {
        res <- data.frame(lambsum = lambsum)
      }
    }else{
      warning("Only input/output orientations are implemented!")
      res <- data.frame(lambsum = lambsum)
    }
    
  }else {
    k <- do.call(rbind, lapply(deamodel$DMU, function(x) x$multiplier_rts))
    dimnames(k)[[2]] <- "k"
    if (deamodel$orientation == "io") {
      switch(deamodel$rts,
             crs = {
               rts <- ifelse(abs(k) > thr , "Error", "Constant")
             },
             vrs = {
               rts <- ifelse(k < -thr, "Decreasing", ifelse(abs(k) < thr, "Constant", "Increasing"))
             },
             nirs = {
               rts <- ifelse(k < -thr, "Decreasing","Constant")
             },
             ndrs = {
               rts <- ifelse(k > thr, "Increasing","Constant")
             },
             stop("General returns to scale not implemented yet!")
      )
    }else{
      switch(deamodel$rts,
             crs = {
               rts <- ifelse(abs(k) > thr, "Error", "Constant")
             },
             vrs = {
               rts <- ifelse(k < -thr, "Increasing", ifelse(abs(k) < thr, "Constant", "Decreasing"))
             },
             nirs = {
               rts <- ifelse(k > thr, "Decreasing", "Constant")
             },
             ndrs = {
               rts <- ifelse(k < -thr, "Increasing", "Constant")
             },
             stop("General returns to scale not implemented yet!")
      )
    }
    
    res <- data.frame(k = k, rts = rts)
    colnames(res) <- c("k", "rts")
   }
  
  return(res)
}