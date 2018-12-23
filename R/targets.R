#' @title Targets
#'   
#' @description Extract the targets of the DMUs from a \code{dea} or \code{dea_fuzzy} solution.
#' 
#' @usage targets(deasol)
#' 
#' @param deasol Object of class \code{dea} or \code{dea_fuzzy} obtained with some of the DEA model functions.
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
#'  data("Coll_Blasco_2006")
#'  data_example <- read_data(Coll_Blasco_2006,
#'                            dmus = 1, 
#'                            ni = 2, 
#'                            no = 2)
#'  result <- model_multiplier(data_example, 
#'                             orientation = "io", 
#'                             rts = "crs")
#'  targets(result)
#'  
#' @export

targets <- function(deasol) {
  
  targetlist <- NULL
  
  if (is.dea(deasol)) {
    
    target_input <- NULL
    target_output <- NULL
    if ("target_input" %in% names(deasol$DMU[[1]])) {
      target_input <- do.call(rbind, lapply(deasol$DMU, function(x)
        x$target_input))
      target_output <- do.call(rbind, lapply(deasol$DMU, function(x)
        x$target_output))
    }
    
    project_input <- NULL
    project_output <- NULL
    if (deasol$modelname %in% c("addsupereff", "sbmsupereff")) {
      project_input <- do.call(rbind, lapply(deasol$DMU, function(x)
        x$project_input))
      project_output <- do.call(rbind, lapply(deasol$DMU, function(x)
        x$project_output))
      targetlist <- list(project_input = project_input,
                         project_output = project_output,
                         target_input = target_input,
                         target_output = target_output)
      
    } else {
      targetlist <- list(target_input = target_input,
                         target_output = target_output)
    }
    
    if(is.null(target_input) && is.null(project_input)) {
      stop("No target/project parameters in this solution!")
    }
    
  } else if (is.dea_fuzzy(deasol)) {
    
    dmunames_eval <- names(deasol$dmu_eval)
    dmunames_ref <- names(deasol$dmu_ref)
    inputnames <- rownames(deasol$data$input$mL)
    outputnames <- rownames(deasol$data$output$mL)
    nde <- length(deasol$dmu_eval)
    ndr <- length(deasol$dmu_ref)
    ni <- length(deasol$data$input$mL[, 1])
    no <- length(deasol$data$output$mL[, 1])
    
    if (grepl("kaoliu", deasol$modelname)) {
      nalpha <- length(deasol$alpha)
      
      target_input.W <- NULL
      target_input.B <- NULL
      target_output.W <- NULL
      target_output.B <- NULL
      if (("target_input" %in% names(deasol$alphacut[[1]]$DMU$Worst[[1]])) &&
          !is.null(deasol$alphacut[[1]]$DMU$Worst[[1]]$target_input)) {
        
        target_input.W <- array(0,
                                dim = c(nde, ni, nalpha),
                                dimnames = list(dmunames_eval, inputnames, names(deasol$alphacut)))
        target_input.B <- target_input.W
        
        for (i in 1:nalpha) {
          target_input.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
            x$target_input))
          target_input.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
            x$target_input))
        }
        
        target_output.W <- array(0,
                                 dim = c(nde, no, nalpha),
                                 dimnames = list(dmunames_eval, outputnames, names(deasol$alphacut)))
        target_output.B <- target_output.W
        
        for (i in 1:nalpha) {
          target_output.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
            x$target_output))
          target_output.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
            x$target_output))
        }
        
        project_input.W <- NULL
        project_output.W <- NULL
        project_input.B <- NULL
        project_output.B <- NULL
        if (grepl("addsupereff", deasol$modelname) || grepl("sbmsupereff", deasol$modelname)) {
          
          project_input.W <- array(0,
                                   dim = c(nde, ni, nalpha),
                                   dimnames = list(dmunames_eval, inputnames, names(deasol$alphacut)))
          project_input.B <- project_input.W
          
          for (i in 1:nalpha) {
            project_input.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
              x$project_input))
            project_input.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
              x$project_input))
          }
          
          project_output.W <- array(0,
                                    dim = c(nde, no, nalpha),
                                    dimnames = list(dmunames_eval, outputnames, names(deasol$alphacut)))
          project_output.B <- project_output.W
          
          for (i in 1:nalpha) {
            project_output.W[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Worst, function(x)
              x$project_output))
            project_output.B[, , i] <- do.call(rbind, lapply(deasol$alphacut[[i]]$DMU$Best, function(x)
              x$project_output))
          }
          targetlist <- list(project_input.W = project_input.W,
                             project_input.B = project_input.B,
                             project_output.W = project_output.W,
                             project_output.B = project_output.B,
                             target_input.W = target_input.W,
                             target_input.B = target_input.B,
                             target_output.W = target_output.W,
                             target_output.B = target_output.B)
          
        } else {
          targetlist <- list(target_input.W = target_input.W,
                             target_input.B = target_input.B,
                             target_output.W = target_output.W,
                             target_output.B = target_output.B)
        }
        
      }
      
      if(is.null(target_input.W) && is.null(project_input.W)) {
        stop("No target/project parameters in this solution!")
      }
      
    } else if (grepl("possibilistic", deasol$modelname)) {
      nh <- length(deasol$h)
      
      if (any(grepl("target", names(deasol$hlevel[[1]]$DMU[[1]])))) {
        
        target_input <- NULL
        if (("target_input" %in% names(deasol$hlevel[[1]]$DMU[[1]])) &&
            !is.null(deasol$hlevel[[1]]$DMU[[1]]$target_input)) {
          
          target_input <- array(0,
                                dim = c(nde, ni, nh),
                                dimnames = list(dmunames_eval, inputnames, names(deasol$hlevel)))
          
          for (i in 1:nh) {
            target_input[, , i] <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x)
              x$target_input))
          }
          
          target_output <- array(0,
                                 dim = c(nde, no, nh),
                                 dimnames = list(dmunames_eval, outputnames, names(deasol$hlevel)))
          
          for (i in 1:nh) {
            target_output[, , i] <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x)
              x$target_output))
          }
          
          targetlist <- list(target_input = target_input,
                             target_output = target_output)
          
        }
        
      } else {
        stop("No target parameters in this solution!")
      }
      
    } else {
      stop("No target parameters in this solution!")
    }
    
  } else {
    
    stop("Input should be a dea or dea_fuzzy class object!")
    
  }
  
  return(targetlist)
  
}
