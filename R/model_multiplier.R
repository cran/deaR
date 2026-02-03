#' @title Multiplier DEA model
#'   
#' @description Solve input-oriented and output-oriented basic DEA  models
#' (multiplicative form) under constant (CCR DEA model), variable (BCC DEA model),
#' non-increasing, non-decreasing or generalized returns to scale. It does not
#' take into account non-controllable, non-discretionary or undesirable inputs/outputs.
#' 
#' @note (1) Very important with the multiplier model: "The optimal weights for
#' an efficient DMU need not be unique" (Cooper, Seiford and Tone, 2007:31).
#' "Usually, the optimal weights for inefficient DMUs are unique, the exception
#' being when the line of the DMU is parallel to one of the boundaries of the
#' feasible region" (Cooper, Seiford and Tone, 2007:32).
#' 
#' (2) The measure of technical input (or output) efficiency obtained by using
#' multiplier DEA models is better the smaller the value of epsilon.
#' 
#' (3) Epsilon is usually set equal to 10^-6. However, if epsilon is not set
#' correctly, the multiplier model can be infeasible (Zhu,2014:49).        
#' 
#' @usage model_multiplier(datadea,
#'                  dmu_eval = NULL,
#'                  dmu_ref = NULL,
#'                  epsilon = 0,
#'                  orientation = c("io", "oo"),
#'                  rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                  L = 1,
#'                  U = 1,
#'                  returnlp = FALSE,
#'                  compute_lambda = TRUE,
#'                  ...)
#' 
#' @param datadea A \code{deadata} object, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param epsilon Numeric, multipliers must be >= \code{epsilon}.
#' @param orientation A string, equal to "io" (input-oriented) or "oo" (output-oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective
#' function and constraints).
#' @param compute_lambda Logical. If it is \code{TRUE}, it computes the dual problem and lambdas.
#' @param ... Ignored, for compatibility issues.
#' 
#' @returns A list of class \code{dea} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{orientation}, \code{rts},
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
#' Charnes, A.; Cooper, W.W. (1962). “Programming with Linear Fractional Functionals”,
#' Naval Research Logistics Quarterly 9, 181-185. \doi{10.1002/nav.3800090303}  
#' 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1978). “Measuring the Efficiency of Decision
#' Making Units”, European Journal of Operational Research 2, 429–444.
#' \doi{10.1016/0377-2217(78)90138-8}  
#' 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1979). “Short Communication: Measuring the
#' Efficiency of Decision Making Units”, European Journal of Operational Research 3, 339.
#' \doi{10.1016/0377-2217(79)90229-7}  
#' 
#' Golany, B.; Roll, Y. (1989). "An Application Procedure for DEA", OMEGA International
#' Journal of Management Science, 17(3), 237-250. \doi{10.1016/0305-0483(89)90029-7}
#' 
#' Seiford, L.M.; Thrall, R.M. (1990). “Recent Developments in DEA. The Mathematical
#' Programming Approach to Frontier Analysis”, Journal of Econometrics 46, 7-38.
#' \doi{10.1016/0304-4076(90)90045-U}
#' 
#' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking.
#' Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York.
#' \doi{10.1007/978-3-319-06647-9}
#' 
#' @examples 
#' # Example 1.
#' # Replication of results in Golany and Roll (1989).
#' data("Golany_Roll_1989")
#' data_example <- make_deadata(datadea = Golany_Roll_1989[1:10, ],
#'                              inputs = 2:4, 
#'                              outputs = 5:6) 
#' result <- model_multiplier(data_example, 
#'                            epsilon = 0, 
#'                            orientation = "io", 
#'                            rts = "crs") 
#' efficiencies(result)
#' multipliers(result)
#' 
#' # Example 2.
#' # Multiplier model with infeasible solutions (See note).
#' data("Fortune500")
#' data_Fortune <- make_deadata(datadea = Fortune500, 
#'                              inputs = 2:4, 
#'                              outputs = 5:6) 
#' result2 <- model_multiplier(data_Fortune, 
#'                            epsilon = 1e-6, 
#'                            orientation = "io", 
#'                            rts = "crs") 
#' # Results for General Motors and Ford Motor are not shown by deaR 
#' # because the solution is infeasible.
#' efficiencies(result2)
#' multipliers(result2)
#' 
#' @seealso \code{\link{model_basic}}, \code{\link{cross_efficiency}}
#' 
#' @import lpSolve
#' 
#' @export
  
model_multiplier <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           epsilon = 0,
           orientation = c("io", "oo"), 
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           returnlp = FALSE,
           compute_lambda = TRUE,
           ...) {
 
  # Checking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
    
  # Checking non-controllable or non-discretionary inputs/outputs
  if ((!is.null(datadea$nc_inputs)) || (!is.null(datadea$nc_outputs))
      || (!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
    warning("This model does not take into account non-controllable or non-discretionary
    feature for inputs/outputs. Instead, you can run model_basic with compute_multipliers = TRUE.")
  }
  
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
  
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
  
  if (rts != "grs") {
    L <- 1
    U <- 1
  } else {
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
  
  if (orientation == "io") {
    input <- datadea$input
    output <- datadea$output
    orient <- 1
  } else {
    input <- -datadea$output
    output <- -datadea$input
    orient <- -1
  }
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  ###########################
  
  obj <- "max"
  
  if (rts == "crs") {
    f.con.rs <- rbind(c(rep(0, ni + no), 1, 0),
                      c(rep(0, ni + no), 0, 1))
    f.dir.rs <- c("=", "=")
    f.rhs.rs <- c(0, 0)
  } else if (rts == "nirs") {
    f.con.rs <- c(rep(0, ni + no), 1, 0)
    f.dir.rs <- "="
    f.rhs.rs <- 0
  }else if (rts == "ndrs") {
    f.con.rs <- c(rep(0, ni + no), 0, 1)
    f.dir.rs <- "="
    f.rhs.rs <- 0
  } else {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  }
  
  if (epsilon > 0) {
    f.con.eps <- cbind(diag(ni + no), matrix(0, nrow = ni + no, ncol = 2))
    f.dir.eps <- rep(">=", ni + no)
    f.rhs.eps <- rep(epsilon, ni + no)
  } else {
    f.con.eps <- NULL
    f.dir.eps <- NULL
    f.rhs.eps <- NULL
  }
    
  # Constraints matrix of 2nd bloc of constraints
  f.con.2 <- cbind(-t(inputref), t(outputref), matrix(1, nrow = ndr, ncol = 1),
                   matrix(-1, nrow = ndr, ncol = 1))
  
  # Directions vector
  f.dir <- c("=", rep("<=", ndr), f.dir.eps, f.dir.rs) # Efficiency is considered free
  #f.dir <- c("<=", rep("<=", ndr), f.dir.eps, f.dir.rs) # Efficiency is considered non-negative
  
  # Right hand side vector
  f.rhs <- c(orient, rep(0, ndr), f.rhs.eps, f.rhs.rs)
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    # Objective function coefficients
    f.obj <- c(rep(0, ni), output[, ii], L, -U)
      
    # Constraints matrix
    f.con.1 <- c(input[, ii], rep(0, no), 0, 0)
    f.con <- rbind(f.con.1, f.con.2, f.con.eps, f.con.rs)
       
    if (returnlp) {
      
      multiplier_input <- rep(0, ni)
      names(multiplier_input) <- inputnames
      multiplier_output <- rep(0, no)
      names(multiplier_output) <- outputnames
      multiplier_rts <- rep(0, 2)
      names(multiplier_rts) <- c("grs_L", "grs_U")
      if (orientation == "io") {
        var = list(multiplier_input = multiplier_input, multiplier_output = multiplier_output,
                   multiplier_rts = multiplier_rts)
      } else {
        var = list(multiplier_output = multiplier_input, multiplier_input = multiplier_output,
                   multiplier_rts = multiplier_rts)
      }
      DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                       const.dir = f.dir, const.rhs = f.rhs, var = var)
      
    } else {
    
      if (compute_lambda) {
        
        res <- lp(obj, f.obj, f.con, f.dir, f.rhs, compute.sens = TRUE)
        
        if (res$status == 0) {
          
          efficiency <- orient * res$objval
          
          lambda <- res$duals[2 : (ndr + 1)]
          names(lambda) <- dmunames[dmu_ref]
          
          target_input <- orient * as.vector(inputref %*% lambda)
          names(target_input) <- inputnames
          target_output <- orient * as.vector(outputref %*% lambda)
          names(target_output) <- outputnames
          
          slack_input <- efficiency * input[, ii] - orient * target_input
          names(slack_input) <- inputnames
          slack_output <- orient * target_output - output[, ii]
          names(slack_output) <- outputnames
          
        } else {
          
          efficiency <- NA
          lambda <- NA
          target_input <- NA
          target_output <- NA
          slack_input <- NA
          slack_output <- NA
          
        }
        
      } else {
        
        res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
        
        if (res$status == 0) {
          efficiency <- orient * res$objval
        } else {
          efficiency <- NA
        }
        
      }
      
      if (res$status == 0) {
        
        multiplier_input <- res$solution[1 : ni]
        names(multiplier_input) <- inputnames
        multiplier_output <- res$solution[(ni + 1) : (ni + no)]
        names(multiplier_output) <- outputnames
        if (rts == "vrs") {
          multiplier_rts <- res$solution[1 + ni + no] - res$solution[2 + ni + no]
          names(multiplier_rts) <- "vrs"
        } else if (rts == "nirs") {
          multiplier_rts <- -res$solution[2 + ni + no]
          names(multiplier_rts) <- "nirs"
        } else if (rts == "ndrs") {
          multiplier_rts <- res$solution[1 + ni + no]
          names(multiplier_rts) <- "ndrs"
        } else if (rts == "grs") {
          multiplier_rts <- c(res$solution[1 + ni + no], -res$solution[2 + ni + no])
          names(multiplier_rts) <- c("grs_L", "grs_U")
        }
        
        if (orientation == "oo") {
          aux <- multiplier_input
          multiplier_input <- multiplier_output
          multiplier_output <- aux
          if (rts != "crs") {
            multiplier_rts <- -multiplier_rts
          }
          # efficiency <- 1 / efficiency # Alternative
        }
        
      } else {
        multiplier_input <- NA
        multiplier_output <- NA
        multiplier_rts <- NA
      }
      
      if (rts == "crs") {
        if (compute_lambda) {
          DMU[[i]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_input, multiplier_output = multiplier_output,
                           lambda = lambda,
                           slack_input = slack_input, slack_output = slack_output,
                           target_input = target_input, target_output = target_output)
        } else {
          DMU[[i]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_input, multiplier_output = multiplier_output)
        }
      } else {
        if (compute_lambda) {
          DMU[[i]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_input, multiplier_output = multiplier_output,
                           multiplier_rts = multiplier_rts,
                           lambda = lambda,
                           slack_input = slack_input, slack_output = slack_output,
                           target_input = target_input, target_output = target_output)
        } else {
          DMU[[i]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_input, multiplier_output = multiplier_output,
                           multiplier_rts = multiplier_rts)
        }
      }

      
      
    }
    
  }
  
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
 
  deaOutput <- list(modelname = "multiplier",
                    orientation = orientation,
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    epsilon = epsilon)
 
  return(structure(deaOutput, class = "dea"))
 
}
