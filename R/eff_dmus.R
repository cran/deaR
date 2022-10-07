#' @title Efficient DMUs.
#' 
#' @description Returns the efficient DMUs evaluated in a \code{dea} class object.
#' 
#' @note If \code{maxslack} is \code{FALSE}, the slacks computed in the first stage
#' are supposed to be the max slacks.
#' 
#' @param deasol An object of class \code{"dea"} obtained by a dea model function.
#' @param tol Numeric. Absolute tolerance for numeric comparisons in efficiency scores.
#' By default, it is 1e-4.
#' 
#' @return A numeric vector containing which DMUs has been evaluated as efficient.
#' This vector is empty if there is not any efficient DMU.
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
#' dataFortune <- read_data(Fortune500,
#'                          ni = 3,
#'                          no = 2)
#' ccrFortune <- model_basic(dataFortune)
#' eff_dmus(ccrFortune)
#' 
#' @export 

eff_dmus <- function(deasol, tol = 1e-4) {

  # Cheking whether deasol is of class "dea" or not...  
  if (!is.dea(deasol)) {
    stop("Data should be of class dea. Run a model first!")
  }
  
  modelname <- deasol$modelname
  orientation <- deasol$orientation
  dmu_eval <- deasol$dmu_eval
  dmu_ref <- deasol$dmu_ref
  nde <- length(dmu_eval)
  iseff <- rep(FALSE, nde)
  
  if (!modelname %in% c("additive", "addsupereff", "basic", "deaps", "fdh_basic", "multiplier",
                        "nonradial", "profit", "rdm", "sbmeff", "sbmsupereff", "supereff_basic")) {
    stop("Model not supported.")
  }
  
  if (modelname %in% c("addsupereff", "sbmsupereff", "supereff_basic")) {
    super <- TRUE
  } else {
    super <- FALSE
  }
  
  # Where is the efficiency score in the deasol$DMU[[1]] list?
  if (modelname == "profit") {
    iesc <- 2
  } else {
    iesc <- 1
  }
  eff <- unlist(lapply(deasol$DMU, function(x) x[[iesc]]))
  i_nona <- which(!is.na(eff))
  
  if (modelname %in% c("additive", "profit", "sbmeff")) { # No revise slacks
    slack_i <- matrix(0, nrow = nde, ncol = 1)
    slack_o <- slack_i
  } else {                                                # Revise slacks
    if ("slack_input" %in% names(deasol$DMU[[1]])) {
      slack_i <- do.call(rbind, lapply(deasol$DMU, function(x) x$slack_input))
    } else {
      slack_i <- matrix(0, nrow = nde, ncol = 1)
    }
    if ("slack_output" %in% names(deasol$DMU[[1]])) {
      slack_o <- do.call(rbind, lapply(deasol$DMU, function(x) x$slack_output))
    } else {
      slack_o <- matrix(0, nrow = nde, ncol = 1)
    }
  }
  
  if (!super) {
    
    # Score for efficiency
    if ((modelname %in% c("additive", "rdm")) || 
        ((modelname %in% c("basic", "fdh_basic")) && (orientation == "dir"))) {
      effscore <- 0 
    } else {
      effscore <- 1
    }
    
    for (i in i_nona) { #### General case ####
      if ((abs(eff[i] - effscore) < tol) && all(abs(slack_i[i, ]) < tol) &&
          all(abs(slack_o[i, ]) < tol)) {
        iseff[i] <- TRUE
      } else if (!dmu_eval[i] %in% dmu_ref) { #### Particular cases ####
        if (modelname == "basic") {
          if ((orientation == "io") && (eff[i] > (1 + tol))) {
            iseff[i] <- TRUE
          } else if ((orientation == "oo") && (eff[i] < (1 - tol))) {
            iseff[i] <- TRUE
          }
        }
      } 
    }
    
  } else { # Super-efficiency
    
    if ((modelname %in% c("addsupereff", "sbmsupereff")) ||
        ((modelname == "supereff_basic") && (orientation == "io"))) {
      
      for (i in i_nona) {
        if (eff[i] > (1 + tol)) { # score > 1
          iseff[i] <- TRUE
        } else if (eff[i] > (1 - tol)) { # score == 1
          if (all(abs(slack_i[i, ]) < tol) && all(abs(slack_o[i, ]) < tol)) {
            iseff[i] <- TRUE
          }
        }
      }
      
    } else { # "supereff_basic" "oo"
      
      for (i in i_nona) {
        if (eff[i] < (1 - tol)) { # score < 1
          iseff[i] <- TRUE
        } else if (eff[i] < (1 + tol)) { # score == 1
          if (all(abs(slack_i[i, ]) < tol) && all(abs(slack_o[i, ]) < tol)) {
            iseff[i] <- TRUE
          }
        }
      }
      
    }
    
  }
  
  res <- dmu_eval[iseff]
  
  return(res)
  
}
