#' @title Non-radial DEA model.
#'   
#' @description Non-radial DEA model allows for non-proportional reductions in each input or augmentations in each output. 
#' 
#' @usage model_nonradial(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             orientation = c("io", "oo"),
#'             rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'             L = 1,
#'             U = 1,
#'             maxslack = TRUE,
#'             weight_slack = 1,
#'             compute_target = TRUE,
#'             returnlp = FALSE,
#'             ...)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' @param orientation A string, equal to "io" (input-oriented) or "oo" (output-oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#'            "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param maxslack Logical. If it is \code{TRUE}, it computes the max slack solution.
#' @param weight_slack If input-oriented, it is a value, vector of length \code{s}, or matrix \code{s} x \code{ne} (where \code{ne} is the length of \code{dmu_eval}) with the weights of the output slacks for the max slack solution.
#'                     If output-oriented, it is a value, vector of length \code{m}, or matrix \code{m} x \code{ne} with the weights of the input slacks for the max slack solution.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets of the max slack solution. 
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective function and constraints) of stage 1.
#' @param ... Ignored, for compatibility issues.
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
#' Banker, R.D.; Morey, R.C. (1986). "Efficiency Analysis for Exogenously Fixed Inputs and Outputs", Operations Research, 34, 80-97. \url{https://doi.org/10.1287/opre.34.4.513}
#' 
#' Färe, R.; Lovell, C.K. (1978). "Measuring the Technical Efficiency of Production", Journal of Economic Theory, 19(1), 150-162. \url{https://doi.org/10.1016/0022-0531(78)90060-1}
#'
#' Wh, J.; Tsai, H.; Zhou, Z. (2011). "Improving Efficiency in International Tourist Hotels in Taipei Using a Non-Radial DEA Model", International Journal of Contemporary Hospitatlity Management, 23(1), 66-83. \url{https://doi.org/10.1108/09596111111101670} 
#'
#' Zhu, J. (1996). “Data Envelopment Analysis with Preference Structure”, The Journal of the Operational Research Society, 47(1), 136. DOI: 10.2307/2584258  
#' 
#' @examples 
#' # Replication of results in Wu, Tsai and Zhou (2011)
#' data("Hotels")
#' data_hotels <- read_data(Hotels, 
#'                          dmus = 1, 
#'                          inputs = 2:5, 
#'                          outputs = 6:8)
#' result <- model_nonradial(data_hotels, 
#'                           orientation = "oo", 
#'                           rts = "vrs")
#' efficiencies(result)
#' 
#' @seealso \code{\link{model_deaps}}, \code{\link{model_sbmeff}}
#'  
#' @import lpSolve
#' 
#' @export
  
model_nonradial <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = c("io", "oo"), 
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           maxslack = TRUE,
           weight_slack = 1,
           compute_target = TRUE,
           returnlp = FALSE,
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
  
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
  # Checking undesirable io and rts
  #if (((!is.null(datadea$ud_inputs)) || (!is.null(datadea$ud_outputs))) && (rts != "vrs")) {
  #  rts <- "vrs"
  #  warning("Returns to scale changed to variable (vrs) because there is data with undesirable inputs/outputs.")
  #}
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
  
  if (rts == "grs") {
    if (L > 1) {
      stop("L must be <= 1.")
    }
    if (U < 1) {
      stop("U must be >= 1.")
    }
  }
  
  dmunames <- datadea$dmunames
  nd <- length(dmunames) # number of dmus
  
  if (is.null(dmu_eval)) {
    dmu_eval <- 1:nd
  } else if (all(dmu_eval %in% (1:nd)) == FALSE) {
    stop("Invalid set of DMUs to be evaluated (dmu_eval).")
  }
  names(dmu_eval) <- dmunames[dmu_eval]
  nde <- length(dmu_eval)
  
  if (is.null(dmu_ref)) {
    dmu_ref <- 1:nd
  } else if (all(dmu_ref %in% (1:nd)) == FALSE) {
    stop("Invalid set of reference DMUs (dmu_ref).")
  }
  names(dmu_ref) <- dmunames[dmu_ref]
  ndr <- length(dmu_ref)
  
  if (orientation == "io") {
    input <- datadea$input
    output <- datadea$output
    nc_inputs <- datadea$nc_inputs
    nc_outputs <- datadea$nc_outputs
    obj <- "min"
    orient <- 1
  } else {
    input <- -datadea$output
    output <- -datadea$input
    nc_inputs <- datadea$nc_outputs
    nc_outputs <- datadea$nc_inputs
    obj <- "max"
    orient <- -1
  }
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  

  nnci <- length(nc_inputs)
  if (nnci == ni) {
    stop("All inputs (io) or outputs (oo) are non-controllable.")
  }
  
  # Checking weights
  if(is.null(weight_slack)){
    weight_slack  <- 1
  }
  if (is.matrix(weight_slack)) {
    if ((nrow(weight_slack) != no) || (ncol(weight_slack) != nde)) {
      stop("Invalid slack weights matrix (number of inputs (io) or outputs (oo) x number of evaluated DMUs).")
    }
  } else if ((length(weight_slack) == 1) || (length(weight_slack) == no)) {
    weight_slack <- matrix(weight_slack, nrow = no, ncol = nde)
  } else {
    stop("Invalid slack weights vector (number of inputs (io) or outputs (oo)).")
  }
  rownames(weight_slack) <- outputnames
  colnames(weight_slack) <- dmunames[dmu_eval]

  target_input <- NULL
  target_output <- NULL
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  ###########################
  
  # Vector de coeficientes de la función objetivo stage 1
  f.obj <- c(rep(1 / (ni - nnci), ni), rep(0, ndr))
  f.obj[nc_inputs] <- 0
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.con2.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- cbind(matrix(0, nrow = 1, ncol = ni), matrix(1, nrow = 1, ncol = ndr))
    f.con2.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = no))
    f.rhs.rs <- 1
    if (rts == "vrs") {
      f.dir.rs <- "="
    } else if (rts == "nirs") {
      f.dir.rs <- "<="
    } else if (rts == "ndrs") {
      f.dir.rs <- ">="
    } else {
      f.con.rs <- rbind(f.con.rs, f.con.rs)
      f.con2.rs <- rbind(f.con2.rs, f.con2.rs)
      f.dir.rs <- c(">=", "<=")
      f.rhs.rs <- c(L, U)
    }
  }
  
  # Matriz técnica del 2º y 3º bloque de restricciones stage 1
  f.con.2 <- cbind(matrix(0, nrow = no, ncol = ni), outputref)
  f.con.3 <- cbind(orient * diag(ni), matrix(0, nrow = ni, ncol = ndr))
  
  if (maxslack && (!returnlp)) {
    
    nnco <- length(nc_outputs) # number of non-controllable outputs
    
    # Matriz técnica stage 2
    f.con2.1 <- cbind(inputref, matrix(0, nrow = ni, ncol = no))
    
    f.con2.2 <- cbind(outputref, -diag(no))
    f.con2.2[nc_outputs, (ndr + 1) : (ndr + no)] <- 0
    
    f.con2.nc <- matrix(0, nrow = nnco, ncol = (ndr + no))
    f.con2.nc[, ndr + nc_outputs] <- diag(nnco)
    
    f.con2 <- rbind(f.con2.1, f.con2.2, f.con2.nc, f.con2.rs)
    
    # Vector de dirección de restricciones stage 2
    f.dir2 <- c(rep("=", ni + no + nnco), f.dir.rs)
    
  }
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    # Matriz técnica stage 1
    f.con.1 <- cbind(-diag(input[, ii], nrow = ni), inputref)
    f.con <- rbind(f.con.1, f.con.2, f.con.3, f.con.rs)
    
    # Vector de dirección de restricciones stage 1
    f.dir <- c(rep("=", ni), rep(">=", no), rep("<=", ni), f.dir.rs)
    f.dir[ni + c(nc_outputs, no + nc_inputs)] <- "="
    
    # Vector de términos independientes stage 1
    f.rhs <- c(rep(0, ni), output[, ii], rep(orient, ni), f.rhs.rs)
    
    if (returnlp) {
      
      efficiency = rep(0, ni)
      names(efficiency) <- inputnames
      lambda <- rep(0, ndr)
      names(lambda) <- dmunames[dmu_ref]
      var <- list(efficiency = efficiency, lambda = lambda)
      DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs,
                       var = var)
      
    } else {
      
      res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
      
      if (res$status == 0) {
        
        mean_eff <- res$objval
        eff <- res$solution[1 : ni]
        names(eff) <- inputnames
        
        if (maxslack) {
          
          # Vector de coeficientes de la función objetivo stage 2
          f.obj2 <- c(rep(0, ndr), weight_slack[, i])
          
          # Vector de términos independientes stage 2
          f.rhs2 <- c(eff * input[, ii], output[, ii], rep(0, nnco), f.rhs.rs)
          
          res <- lp("max", f.obj2, f.con2, f.dir2, f.rhs2)$solution
          
          lambda <- res[1 : ndr]
          names(lambda) <- dmunames[dmu_ref]
          
          slack_output <- res[(ndr + 1) : (ndr + no)]
          names(slack_output) <- outputnames
          
          if (compute_target) {
            target_input <- orient * as.vector(inputref %*% lambda)
            target_output <- orient * as.vector(outputref %*% lambda)
            #target_input <- orient * eff * input[, ii]
            names(target_input) <- inputnames
            #target_output <- orient * (output[, ii] + slack_output)
            names(target_output) <- outputnames
          }
          
        } else {
          
          lambda <- res$solution[(ni + 1) : (ni + ndr)]
          names(lambda) <- dmunames[dmu_ref]
          
          target_input <- orient * as.vector(inputref %*% lambda)
          #target_input <- orient * eff * input[, ii]
          names(target_input) <- inputnames
          target_output <- orient * as.vector(outputref %*% lambda)
          names(target_output) <- outputnames
          
          slack_output <- orient * target_output - output[, ii]
          names(slack_output) <- outputnames
          
        }
        
      } else {
        
        mean_eff <- NA
        eff <- NA
        lambda <- NA
        slack_output <- NA
        if (compute_target) {
          target_input <- NA
          target_output <- NA
        }
        
      }
      
      if (orientation == "io") {
        DMU[[i]] <- list(mean_efficiency = mean_eff,
                         efficiency = eff,
                         lambda = lambda,
                         slack_output = slack_output,
                         target_input = target_input, target_output = target_output)
      } else {
        DMU[[i]] <- list(mean_efficiency = mean_eff,
                         efficiency = eff,
                         lambda = lambda,
                         slack_input = slack_output,
                         target_input = target_output, target_output = target_input)
      }
      
    }
    
  }
 
  deaOutput <- list(modelname = "nonradial",
                    orientation = orientation,
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref)
 
  return(structure(deaOutput, class = "dea"))
 
}
