#' @title Efficiencies
#'   
#' @description Extract the scores (optimal objective values) of the evaluated DMUs
#' from a DEA fuzzy solution. Note that these scores may not always be interpreted
#' as efficiencies.
#' 
#' @param x Object of class dea or dea_fuzzy obtained with some of the dea model functions.
#' @param ... Other options (for compatibiliy)
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
#' # Replication results model DEA1 in Tomkins and Green (1988)
#' data("Departments")
#' # Calculate Total income
#' Departments$Total_income <- Departments[, 5] + Departments[, 6] + Departments[, 7] 
#' data_DEA1 <- read_data(Departments,
#'                        inputs = 9,
#'                        outputs = c(2, 3, 4, 12))
#' result <- model_basic(data_DEA1, 
#'                       orientation = "io",
#'                       rts = "crs")
#' efficiencies(result) # Table 3 (p.156) 
#'  
#' @export

efficiencies.dea_fuzzy <-
  function(x, ...) {
    deasol <- x
    dmunames_eval <- names(deasol$dmu_eval)
    nde <- length(deasol$dmu_eval)
    
    if (grepl("kaoliu", deasol$modelname)) {
      nalpha <- length(deasol$alpha)
      
      if ("efficiency" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) {
        
        neff <- length(deasol$alphacut[[1]]$DMU$Worst[[1]]$efficiency)
        
         if (neff == 1) {

          eff.W <- matrix(0, nrow = nde, ncol = nalpha)
          rownames(eff.W) <- dmunames_eval
          colnames(eff.W) <- names(deasol$alphacut)
          eff.B <- eff.W

          for (j in 1:nalpha) {
            eff.W[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Worst, function(x)
              x$efficiency))
            eff.B[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Best, function(x)
              x$efficiency))
          }
          
        } else {
          
          eff.W <- array(0,
                         dim = c(nde, neff + 1, nalpha),
                         dimnames = list(dmunames_eval,
                                         c(names(deasol$alphacut[[1]]$DMU$Worst[[1]]$efficiency), "mean_efficiency"),
                                         names(deasol$alphacut)))
          eff.B <- eff.W
          
          for (k in 1:nalpha) {
            eff.W[, , k]  <- cbind(
              do.call(rbind, lapply(deasol$alphacut[[k]]$DMU$Worst, function(x)
              x$efficiency)),
              unlist(lapply(deasol$alphacut[[k]]$DMU$Worst, function(x)
                x$mean_efficiency))
            )
            eff.B[, , k]  <- cbind(
              do.call(rbind, lapply(deasol$alphacut[[k]]$DMU$Best, function(x)
                x$efficiency)),
              unlist(lapply(deasol$alphacut[[k]]$DMU$Best, function(x)
                x$mean_efficiency))
            )
          }
          
        }
        
      } else if ("beta" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) {
        
        eff.W <- matrix(0, nrow = nde, ncol = nalpha)
        rownames(eff.W) <- dmunames_eval
        colnames(eff.W) <- names(deasol$alphacut)
        eff.B <- eff.W
        
        for (j in 1:nalpha) {
          eff.W[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Worst, function(x)
            x$beta))
          eff.B[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Best, function(x)
            x$beta))
        }
        
      } else if ("delta" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) {
        
        eff.W <- matrix(0, nrow = nde, ncol = nalpha)
        rownames(eff.W) <- dmunames_eval
        colnames(eff.W) <- names(deasol$alphacut)
        eff.B <- eff.W
        
        for (j in 1:nalpha) {
          eff.W[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Worst, function(x)
            x$delta))
          eff.B[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Best, function(x)
            x$delta))
        }
        
      } else if ("objval" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) {
        
        eff.W <- matrix(0, nrow = nde, ncol = nalpha)
        rownames(eff.W) <- dmunames_eval
        colnames(eff.W) <- names(deasol$alphacut)
        eff.B <- eff.W
        
        for (j in 1:nalpha) {
          eff.W[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Worst, function(x)
            x$objval))
          eff.B[, j] <- unlist(lapply(deasol$alphacut[[j]]$DMU$Best, function(x)
            x$objval))
        }
        
      } else {
        stop("No efficiency/beta/delta/objval parameters in this solution!")
      }
      
      return(list(Worst = round(eff.W, 5), Best = round(eff.B, 5)))
      
    } else if (grepl("possibilistic", deasol$modelname)) {
      nh <- length(deasol$h)
      
      if ("efficiency" %in% names(deasol$hlevel[[1]]$DMU[[1]])) {
        
        neff <- length(deasol$hlevel[[1]]$DMU[[1]]$efficiency)
        
        if (neff == 1) {
          
          eff <- matrix(0, nrow = nde, ncol = nh)
          rownames(eff) <- dmunames_eval
          colnames(eff) <- names(deasol$hlevel)
          
          for (j in 1:nh) {
            eff[, j] <- unlist(lapply(deasol$hlevel[[j]]$DMU, function(x)
              x$efficiency))
          }
          
        } else {
          
          eff <- array(0,
                       dim = c(nde, neff + 1, nh),
                       dimnames = list(dmunames_eval,
                                       c(names(deasol$hlevel[[1]]$DMU[[1]]$efficiency), "mean_efficiency"),
                                       names(deasol$hlevel)))
          
          for (k in 1:nh) {
            eff[, , k]  <- cbind(
              do.call(rbind, lapply(deasol$hlevel[[k]]$DMU, function(x)
                x$efficiency)),
              unlist(lapply(deasol$hlevel[[k]]$DMU, function(x)
                x$mean_efficiency))
            )
          }
          
        }
        
      } else if ("beta" %in% names(deasol$hlevel[[1]]$DMU[[1]])) {
        
        eff <- matrix(0, nrow = nde, ncol = nh)
        rownames(eff) <- dmunames_eval
        colnames(eff) <- names(deasol$hlevel)
        
        for (j in 1:nh) {
          eff[, j] <- unlist(lapply(deasol$hlevel[[j]]$DMU, function(x)
            x$beta))
        }
      
      } else {
        stop("No efficiency/beta parameters in this solution!")
      }
      
      return(round(eff, 5))
      
    } else if (grepl("guotanaka", deasol$modelname)) {
      nh <- length(deasol$h)
      
      eff <- array(0,
                   dim = c(nde, 3, nh),
                   dimnames = list(dmunames_eval,
                                   names(deasol$hlevel[[1]]$DMU[[1]]$efficiency),
                                   names(deasol$hlevel)))
      
      for (k in 1:nh) {
        eff[, , k]  <- do.call(rbind, lapply(deasol$hlevel[[k]]$DMU, function(x)
            x$efficiency))
      }
      
      return(round(eff, 5))
      
    }
  
}