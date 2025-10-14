#' @title Linear Generalized Oriented DEA model.
#'   
#' @description It solves linear generalized oriented DEA models (see Bolós et
#' al. 2026). By default, models are solved in a two-stage process (slacks are
#' maximized).  
#' 
#' @usage model_lgo(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             d_input = 1,
#'             d_output = 1,
#'             rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'             L = 1,
#'             U = 1,
#'             maxslack = TRUE,
#'             weight_slack_i = 1,
#'             weight_slack_o = 1,
#'             returnlp = FALSE,
#'             ...)
#' 
#' @param datadea A \code{deadata} object with \code{n} DMUs, \code{m} inputs and \code{s}
#' outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation
#' reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param d_input A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the input orientation parameters.
#' If \code{d_input} == 1 (default) and \code{d_output} == 0, it is equivalent
#' to input oriented.
#' @param d_output A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the output orientation parameters.
#' If \code{d_input} == 0 and \code{d_output} == 1 (default), it is equivalent
#' to output oriented.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param maxslack Logical. If it is \code{TRUE}, it computes the max slack solution.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the input slacks
#' for the max slack solution.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the output
#' slacks for the max slack solution.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems
#' (objective function and constraints) of stage 1.
#' @param ... Ignored, for compatibility issues.
#' 
#' @returns A list of class \code{dea} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{orientation_param}, \code{rts},
#'  \code{dmu_eval} and \code{dmu_ref}.
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
#' 
#' Bolós, V.J.; Benítez, R.; Coll-Serrano, V (2026). "A new family of models
#' with generalized orientation in data envelopment analysis". International
#' Transactions in Operational Research. \doi{10.1111/itor.70063}
#' 
#' Chambers, R.G.; Chung, Y.; Färe, R. (1996). "Benefit and Distance Functions",
#' Journal of Economic Theory, 70(2), 407-419. 
#' 
#' Chambers, R.G.; Chung, Y.; Färe, R. (1998). "Profit Directional Distance
#' Functions and Nerlovian Efficiency", Journal of Optimization Theory and
#' Applications, 95, 351-354.
#' 
#' @examples 
#' 
#' data("PFT1981") 
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ] 
#' PFT <- make_deadata(PFT, 
#'                     inputs = 2:6, 
#'                     outputs = 7:9 )
#' eval_pft <- model_lgo(PFT)
#' efficiencies(eval_pft)
#'  
#' @seealso \code{\link{model_basic}}, \code{\link{model_dir}}, \code{\link{model_qgo}}
#' 
#' @import lpSolve
#' 
#' @export
  
model_lgo <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           d_input = 1,
           d_output = 1,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           maxslack = TRUE,
           weight_slack_i = 1,
           weight_slack_o = 1,
           returnlp = FALSE,
           ...) {
 
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
  
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
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
  } else if (!all(dmu_eval %in% (1:nd))) {
    stop("Invalid set of DMUs to be evaluated (dmu_eval).")
  }
  names(dmu_eval) <- dmunames[dmu_eval]
  nde <- length(dmu_eval)
  
  if (is.null(dmu_ref)) {
    dmu_ref <- 1:nd
  } else if (!all(dmu_ref %in% (1:nd))) {
    stop("Invalid set of reference DMUs (dmu_ref).")
  }
  names(dmu_ref) <- dmunames[dmu_ref]
  ndr <- length(dmu_ref)
  
  input <- datadea$input
  output <- datadea$output
  nc_inputs <- datadea$nc_inputs
  nc_outputs <- datadea$nc_outputs
  nd_inputs <- datadea$nd_inputs
  nd_outputs <- datadea$nd_outputs
  ud_inputs <- datadea$ud_inputs
  ud_outputs <- datadea$ud_outputs
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  obj <- "max"
  
  if (is.matrix(d_input)) {
    if ((nrow(d_input) != ni) || (ncol(d_input) != nde)) {
      stop("Invalid input orientation matrix (number of inputs x number of evaluated DMUs).")
    }
  } else if ((length(d_input) == 1) || (length(d_input) == ni)) {
    d_input <- matrix(d_input, nrow = ni, ncol = nde)
  } else {
    stop("Invalid input orientation vector.")
  }
  if (any(d_input < 0)) {
    stop("Input orientation parameters must be non-negative.")
  }
  rownames(d_input) <- inputnames
  colnames(d_input) <- dmunames[dmu_eval]
  
  if (is.matrix(d_output)) {
    if ((nrow(d_output) != no) || (ncol(d_output) != nde)) {
      stop("Invalid output orientation matrix (number of outputs x number of evaluated DMUs).")
    }
  } else if ((length(d_output) == 1) || (length(d_output) == no)) {
    d_output <- matrix(d_output, nrow = no, ncol = nde)
  } else {
    stop("Invalid output orientation vector.")
  }
  if (any(d_output < 0)) {
    stop("Output orientation parameters must be non-negative.")
  }
  rownames(d_output) <- outputnames
  colnames(d_output) <- dmunames[dmu_eval]
  
  dir_input <- d_input * matrix(input[, dmu_eval], nrow = ni) # input of DMUs in dmu_eval
  rownames(dir_input) <- inputnames
  colnames(dir_input) <- dmunames[dmu_eval]
  dir_output <- d_output * matrix(output[, dmu_eval], nrow = no) # output of DMUs in dmu_eval
  rownames(dir_output) <- outputnames
  colnames(dir_output) <- dmunames[dmu_eval]
  
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  ncd_inputs <- c(nc_inputs, nd_inputs)
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  ###########################
  
  # Objective function coefficients stage 1
  f.obj <- c(1, rep(0, ndr))
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.con2.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- cbind(0, matrix(1, nrow = 1, ncol = ndr))
    f.con2.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
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
  
  # Directions vector stage 1
  f.dir <- c(rep("<=", ni), rep(">=", no), f.dir.rs)
  f.dir[c(nc_inputs, ni + nc_outputs)] <- "="
  
  if (maxslack && (!returnlp)) {
    
    # Checking weights
    if (is.matrix(weight_slack_i)) {
      if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
        stop("Invalid weight input matrix (number of inputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
      weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
    } else {
      stop("Invalid weight input vector (number of inputs).")
    }
    rownames(weight_slack_i) <- inputnames
    colnames(weight_slack_i) <- dmunames[dmu_eval]
    weight_slack_i[nd_inputs, ] <- 0 # Non-discretionary io not taken into account for maxslack solution
    
    if (is.matrix(weight_slack_o)) {
      if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
        stop("Invalid weight output matrix (number of outputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
      weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
    } else {
      stop("Invalid weight output vector (number of outputs).")
    }
    rownames(weight_slack_o) <- outputnames
    colnames(weight_slack_o) <- dmunames[dmu_eval]
    weight_slack_o[nd_outputs, ] <- 0 # Non-discretionary io not taken into account for maxslack solution
    
    nnci <- length(nc_inputs) # number of non-controllable inputs
    nnco <- length(nc_outputs) # number of non-controllable outputs
    
    # Constraints matrix stage 2
    f.con2.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = no))
    f.con2.1[nc_inputs, (ndr + 1) : (ndr + ni)] <- 0
    
    f.con2.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no))
    f.con2.2[nc_outputs, (ndr + ni + 1) : (ndr + ni + no)] <- 0
    
    f.con2.nc <- matrix(0, nrow = (nnci + nnco), ncol = (ndr + ni + no))
    f.con2.nc[, ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
    
    f.con2 <- rbind(f.con2.1, f.con2.2, f.con2.nc, f.con2.rs)
    
    # Directions vector stage 2
    f.dir2 <- c(rep("=", ni + no + nnci + nnco), f.dir.rs)
    
    ### Slacks for undesirable i/o are 0 ###
    nudi <- length(ud_inputs)
    nudo <- length(ud_outputs)
    if ((nudi + nudo) > 0) {
      f.con2.s0.1 <- cbind(matrix(0, nrow = nudi, ncol = ndr),
                           matrix(diag(ni)[ud_inputs, ], nrow = nudi, ncol = ni),
                           matrix(0, nrow = nudi, ncol = no))
      f.con2.s0.2 <- cbind(matrix(0, nrow = nudo, ncol = ndr + ni),
                           matrix(diag(no)[ud_outputs, ], nrow = nudo, ncol = no))
      f.con2 <- rbind(f.con2, f.con2.s0.1, f.con2.s0.2)
      f.dir2 <- c(f.dir2, rep("=", nudi + nudo))
    }
    
  }
  
  ncd_outputs <- c(nc_outputs, nd_outputs)
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    # Constraints matrix stage 1
    f.con.1 <- cbind(dir_input[, i], inputref)
    f.con.1[ud_inputs, 1] <- -dir_input[ud_inputs, i]
    f.con.1[ncd_inputs, 1] <- 0
    f.con.2 <- cbind(-dir_output[, i], outputref)
    f.con.2[ud_outputs, 1] <- dir_output[ud_outputs, i]
    f.con.2[ncd_outputs, 1] <- 0
    f.con <- rbind(f.con.1, f.con.2, f.con.rs)
    
    # Directions vector
    f.dir[c(ud_inputs, ni + ud_outputs)] <- "="
    
    # Right hand side vector stage 1
    f.rhs <- c(input[, ii], output[, ii], f.rhs.rs)
    
    if (returnlp) {
      
      lambda <- rep(0, ndr)
      names(lambda) <- dmunames[dmu_ref]
      var <- list(efficiency = 0, lambda = lambda)
      DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                       const.dir = f.dir, const.rhs = f.rhs, var = var)
      
    } else {
      
      res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
      
      if (res$status == 0) {
        
        res <- res$solution
        
        beta <- res[1]
        
        if (maxslack) {
          
          # Objective function coefficients stage 2
          f.obj2 <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i])
          
          # Right hand side vector stage 2
          f.rhs2 <- c(input[, ii] - beta * dir_input[, i], output[, ii] + beta * dir_output[, i],
                      rep(0, nnci + nnco), f.rhs.rs, rep(0, nudi + nudo))
          f.rhs2[ud_inputs] <- input[ud_inputs, ii] + beta * dir_input[ud_inputs, i]
          f.rhs2[ni + ud_outputs] <- output[ud_outputs, ii] - beta * dir_output[ud_outputs, i]
          f.rhs2[ncd_inputs] <- input[ncd_inputs, ii]
          f.rhs2[ni + ncd_outputs] <- output[ncd_outputs, ii]
          
          res <- lp("max", f.obj2, f.con2, f.dir2, f.rhs2)$solution
          
          lambda <- res[1 : ndr]
          names(lambda) <- dmunames[dmu_ref]
          
          slack_input <- res[(ndr + 1) : (ndr + ni)]
          names(slack_input) <- inputnames
          slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
          names(slack_output) <- outputnames
          
          effproj_input <- as.vector(inputref %*% lambda)
          effproj_output <- as.vector(outputref %*% lambda)
          names(effproj_input) <- inputnames
          names(effproj_output) <- outputnames
          
        } else {
          
          lambda <- res[2 : (ndr + 1)]
          names(lambda) <- dmunames[dmu_ref]
          
          effproj_input <- as.vector(inputref %*% lambda)
          names(effproj_input) <- inputnames
          effproj_output <- as.vector(outputref %*% lambda)
          names(effproj_output) <- outputnames
          
          slack_input <- input[, ii] - beta * dir_input[, i] - effproj_input
          slack_input[ud_inputs] <- input[ud_inputs, ii] + beta * dir_input[ud_inputs, i] -
            effproj_input[ud_inputs]
          slack_input[ncd_inputs] <- input[ncd_inputs, ii] - effproj_input[ncd_inputs]
          names(slack_input) <- inputnames
          slack_output <- effproj_output - output[, ii] - beta * dir_output[, i]
          slack_output[ud_outputs] <- effproj_output[ud_outputs] - output[ud_outputs, ii] +
            beta * dir_output[ud_outputs, i]
          slack_output[ncd_outputs] <- effproj_output[ncd_outputs] - output[ncd_outputs, ii]
          names(slack_output) <- outputnames
          
        }
        
        target_input <- effproj_input + slack_input
        target_output <- effproj_output - slack_output
        names(target_input) <- inputnames
        names(target_output) <- outputnames
        
        rho <- (1 - sum(dir_input[, i] / input[, ii]) * beta / ni) /
          (1 + beta * sum(dir_output[, i] / output[, ii]) / no)
        names(rho) <- "rho"
        
      } else {
        
        rho <- NA
        beta <- NA
        lambda <- NA
        target_input <- NA
        target_output <- NA
        slack_input <- NA
        slack_output <- NA
        effproj_input <- NA
        effproj_output <- NA
        
      }
      
      DMU[[i]] <- list(efficiency = rho,
                       beta = beta,
                       lambda = lambda,
                       target_input = target_input, target_output = target_output,
                       slack_input = slack_input, slack_output = slack_output,
                       effproj_input = effproj_input, effproj_output = effproj_output
      )
      
    }
    
  }
  
  orientation_param <- list(
    d_input = d_input,
    d_output = d_output)
  
  # Checking if a DMU is in its own reference set (when rts = "grs")
  if (rts == "grs") {
    eps <- 1e-6
    for (i in 1:nde) {
      j <- which(dmu_ref == dmu_eval[i])
      if (length(j) == 1) {
        kk <- DMU[[i]]$lambda[j]
        kk2 <- sum(DMU[[i]]$lambda[-j])
        if ((kk > eps) && (kk2 > eps)) {
          warning(paste("Under generalized returns to scale,", dmunames[dmu_eval[i]],
                        "appears in its own reference set."))
        }
      }
    }
  }
  
  deaOutput <- list(modelname = "lgo",
                   orientation_param = orientation_param,
                   rts = rts,
                   L = L,
                   U = U,
                   DMU = DMU,
                   data = datadea,
                   dmu_eval = dmu_eval,
                   dmu_ref = dmu_ref,
                   maxslack = maxslack,
                   weight_slack_i = weight_slack_i,
                   weight_slack_o = weight_slack_o)
 
  return(structure(deaOutput, class = "dea"))
}
