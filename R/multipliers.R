#' @title Multipliers
#'   
#' @description Extract the multipliers of the DMUs from a \code{dea} or
#' \code{dea_fuzzy} solution.
#' 
#' @usage multipliers(deasol)
#' 
#' @param deasol Object of class \code{dea} or \code{dea_fuzzy} obtained with
#' some of the DEA model functions.
#' 
#' @returns A list with the multipliers of the evaluated DMUs.
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
#' multipliers(result)
#'  
#' @export

multipliers <- 
  function(deasol) {
    
    if (is.dea(deasol)) {
      
      if (any(grepl("multiplier", names(deasol$DMU[[1]])))) {
        
        multiplier_input <- NULL
        if ("multiplier_input" %in% names(deasol$DMU[[1]])) {
          multiplier_input <- do.call(rbind, lapply(deasol$DMU, function(x)
            x$multiplier_input))
        }
        
        multiplier_output <- NULL
        if ("multiplier_output" %in% names(deasol$DMU[[1]])) {
          multiplier_output <- do.call(rbind, lapply(deasol$DMU, function(x)
            x$multiplier_output))
        }
        
        multiplier_rts <- NULL
        if ("multiplier_rts" %in% names(deasol$DMU[[1]])) {
          multiplier_rts <- unlist(lapply(deasol$DMU, function(x)
            x$multiplier_rts))
          
          multlist <- list(multiplier_input = round(multiplier_input,6), 
                           multiplier_output = round(multiplier_output,6), 
                           multiplier_rts = multiplier_rts)
        } else {
          
          multlist <- list(multiplier_input = round(multiplier_input,6), 
                           multiplier_output = round(multiplier_output,6))
          
        }
        
      } else {
        stop("No multiplier parameters in this solution!")
      }
      
    }  else if (is.dea_fuzzy(deasol)) {
      
      dmunames_eval <- deasol$data$dmunames[deasol$dmu_eval]
      dmunames_ref <- deasol$data$dmunames[deasol$dmu_ref]
      inputnames <- rownames(deasol$data$input$mL)
      outputnames <- rownames(deasol$data$output$mL)
      nde <- length(deasol$dmu_eval)
      #ndr <- length(deasol$dmu_ref)
      ni <- length(deasol$data$input$mL[, 1])
      no <- length(deasol$data$output$mL[, 1])
      
      if (grepl("kaoliu", deasol$modelname)) {
        nalpha <- length(deasol$alpha)
        
        multiplier_input.W <- NULL
        multiplier_input.B <- NULL
        if (("multiplier_input" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) &&
            !is.null(deasol$alphacut[[1]]$DMU$Worst[[1]]$multiplier_input)) {
          
          multiplier_input.W <- array(0,
                                      dim = c(nde, ni, nalpha),
                                      dimnames = list(dmunames_eval, inputnames, names(deasol$alphacut)))
          multiplier_input.B <- multiplier_input.W
          
          for (i in 1:nalpha) {
            multiplier_input.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
              x$multiplier_input))
            multiplier_input.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
              x$multiplier_input))
          }
        }
        
        multiplier_output.W <- NULL
        multiplier_output.B <- NULL
        if (("multiplier_output" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) &&
            !is.null(deasol$alphacut[[1]]$DMU$Worst[[1]]$multiplier_output)) {
          
          multiplier_output.W <- array(0,
                                       dim = c(nde, no, nalpha),
                                       dimnames = list(dmunames_eval, outputnames, names(deasol$alphacut)))
          multiplier_output.B <- multiplier_output.W
          
          for (i in 1:nalpha) {
            multiplier_output.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
              x$multiplier_output))
            multiplier_output.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
              x$multiplier_output))
          }
          
        }
        
        multiplier_rts.W <- NULL
        multiplier_rts.B <- NULL
        if (("multiplier_rts" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) &&
            !is.null(deasol$alphacut[[1]]$DMU$Worst[[1]]$multiplier_rts)) {
          
          multiplier_rts.W <- array(0,
                                    dim = c(nde, nalpha),
                                    dimnames = list(dmunames_eval, names(deasol$alphacut)))
          multiplier_rts.B <- multiplier_rts.W
          
          for (i in 1:nalpha) {
            multiplier_rts.W[, i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
              x$multiplier_rts))
            multiplier_rts.B[, i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
              x$multiplier_rts))
          }
          
          multlist <- list(multiplier_input.W = multiplier_input.W,
                           multiplier_input.B = multiplier_input.B,
                           multiplier_output.W = multiplier_output.W,
                           multiplier_output.B = multiplier_output.B,
                           multiplier_rts.W = multiplier_rts.W,
                           multiplier_rts.B = multiplier_rts.B)
          
        } else {
          
          multlist <- list(multiplier_input.W = multiplier_input.W,
                           multiplier_input.B = multiplier_input.B,
                           multiplier_output.W = multiplier_output.W,
                           multiplier_output.B = multiplier_output.B)
          
        }
        
        if(is.null(multiplier_input.W) && is.null(multiplier_output.W)) {
          stop("No multiplier parameters in this solution!")
        }
        
      } else if (grepl("guotanaka", deasol$modelname)) {
        nh <- length(deasol$h)
        
        multiplier_input <- array(0,
                                  dim = c(nde, ni, nh),
                                  dimnames = list(dmunames_eval, inputnames, names(deasol$hlevel)))
        for (i in 1:nh) {
          multiplier_input[, , i] <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x)
            x$multiplier_input))
        }
        
        multiplier_output <- array(0,
                                   dim = c(nde, no, nh),
                                   dimnames = list(dmunames_eval, outputnames, names(deasol$hlevel)))
        for (i in 1:nh) {
          multiplier_output[, , i] <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x)
            x$multiplier_output))
        }
        
        multlist <- list(multiplier_input = multiplier_input,
                         multiplier_output = multiplier_output)
        
      } else {
        stop("No multiplier parameters in this solution!")
      }
      
    } else {
      
      stop("Input should be a dea or dea_fuzzy class object!")
      
    }
    
    return(multlist)
    
  }