#' @title Additive-min DEA model.
#'   
#' @description Solve the weighted version of the additive-min (mADD) model of
#' Aparicio et. al (2007) with different returns to scale. For non constant returns
#' to scale, a modification given by Zhu et al. (2018) is done.
#' 
#' @note In this model, the efficiency score is the sum of the slacks. Therefore,
#' a DMU is efficient when the objective value (\code{objval}) is zero.
#' 
#' @usage model_addmin(datadea,
#'                dmu_eval = NULL,
#'                dmu_ref = NULL,
#'                orientation = NULL,
#'                weight_slack_i = 1,
#'                weight_slack_o = 1,
#'                rts = c("crs", "vrs", "nirs", "ndrs"),
#'                method = c("mf", "milp"),
#'                extreff = NULL,
#'                M_d = NULL,
#'                M_lambda = 1e3,
#'                maxfr = NULL,
#'                tol = 1e-6,
#'                silent = TRUE,
#'                compute_target = TRUE,
#'                check_target = FALSE,
#'                returnlp = FALSE,
#'                ...)
#' 
#' @param datadea A \code{deadata} object with \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation This parameter is either \code{NULL} (default) or a string, equal to
#' "io" (input-oriented) or "oo" (output-oriented). It is used to modify the weight slacks. 
#' If input-oriented, \code{weight_slack_o} are taken 0.
#' If output-oriented, \code{weight_slack_i} are taken 0.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval})
#' with the weights of the input slacks. If 0, output-oriented.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval})
#' with the weights of the output slacks. If 0, input-oriented.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing) or "ndrs" (non-decreasing). Under non-increasing
#' or non-decreasing returns to scale, you may set \code{check_target = TRUE} because
#' methods are not reliable. Generalized returns to scale are not available.
#' @param method A string with the method: "mf" (default) for maximal friends, or "milp"
#' for the mixed integer linear program of Aparicio et al. (2007). MILP method is
#' faster but very problematic numerically.
#' @param extreff A vector with the extreme efficient DMUs for "milp" method, as it
#' is returned by function \code{extreme_efficient}.  If \code{NULL} (default)
#' this vector is computed internally.
#' @param M_d Numeric, a big positive quantity for "milp" method. It is an upper
#' bound for auxiliary variables named "d" in Aparicio (2007). If \code{NULL}
#' (default), it is estimated automatically. A very big value can produce catastrophic
#' cancellations. If the results are not correct or the solver hangs, try to change its value.
#' @param M_lambda Numeric, a big positive quantity for "milp" method. It is an upper
#' bound for lambda variables. A very big value can produce catastrophic cancellations.
#' If the results are not correct or the solver hangs, try to change its value (1e3 by default).
#' @param maxfr A list with the maximal friends sets for "mf" method, as it is returned by function
#' \code{maximal_friends}. If \code{NULL} (default) this list is computed internally.
#' @param tol Numeric, a tolerance margin for checking efficiency in \code{extreme_efficient}
#' or \code{maximal_friends} functions, and for checking targets.
#' @param silent Logical. If \code{FALSE}, it prints all the messages from
#' function \code{maximal_friends}.
#' @param compute_target Logical. If it is \code{TRUE} (default), it computes targets. 
#' @param check_target Logical. If it is \code{TRUE}, it checks the efficiency of targets.
#' If a target is not efficient, the method has failed.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems
#' (objective function and constraints).
#' @param ... For compatibility issues.
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
#' Aparicio, J.; Ruiz, J.L.; Sirvent, I. (2007) "Closest targets and minimum
#' distance to the Pareto-efficient frontier in DEA", Journal of Productivity
#' Analysis, 28, 209-218. \doi{10.1007/s11123-007-0039-5}
#' 
#' Zhu, Q.; Wu, J.; Ji, X.; Li, F. (2018) "A simple MILP to determine closest
#' targets in non-oriented DEA model satisfying strong monotonicity", Omega, 79,
#' 1-8. \doi{10.1016/j.omega.2017.07.003}
#' 
#' @examples
#' # Example 1.
#' data("Airlines")
#' datadea <- make_deadata(Airlines,
#'                         inputs = 4:7,
#'                         outputs = 2:3)
#' result <- model_addmin(datadea = datadea,
#'                        method = "milp")
#' targets(result)
#' 
#' \dontrun{
#' # Example 2. Directional model with Additive-min model in second stage 
#' data("Airlines")
#' datadea <- make_deadata(Airlines,
#'                         inputs = 4:7,
#'                         outputs = 2:3)
#' resdir <- model_basic(datadea = datadea,
#'                       orientation = "dir",
#'                       maxslack = FALSE)
#' proj_input <- targets(resdir)[[1]] + slacks(resdir)[[1]]
#' proj_output <- targets(resdir)[[2]] - slacks(resdir)[[2]]
#' nd <- ncol(datadea$dmunames) # Number of DMUs
#' maxfr <- maximal_friends(datadea = datadea)
#' for (i in 1:nd) {
#'   datadea2 <- datadea
#'   datadea2$input[, i] <- proj_input[i, ]
#'   datadea2$output[, i] <- proj_output[i, ]
#'   DMUaux <- model_addmin(datadea = datadea2,
#'                          method = "mf",
#'                          maxfr = maxfr,
#'                          dmu_eval = i)$DMU[[1]]
#'   resdir$DMU[[i]]$slack_input <- DMUaux$slack_input
#'   resdir$DMU[[i]]$slack_output <- DMUaux$slack_output
#'   resdir$DMU[[i]]$target_input <- DMUaux$target_input
#'   resdir$DMU[[i]]$target_output <- DMUaux$target_output
#' }
#' targets(resdir)
#' }
#' 
#' @seealso \code{\link{model_additive}}, \code{\link{extreme_efficient}},
#' \code{\link{maximal_friends}}
#' 
#' @import lpSolve
# Rsymphony
#' 
#' @export

model_addmin <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = NULL,
           weight_slack_i = 1,
           weight_slack_o = 1,
           rts = c("crs", "vrs", "nirs", "ndrs"),
           method = c("mf", "milp"),
           extreff = NULL,
           M_d = NULL,
           M_lambda = 1e3,
           maxfr = NULL,
           tol = 1e-6,
           silent = TRUE,
           compute_target = TRUE,
           check_target = FALSE,
           returnlp = FALSE,
           # Rsym = FALSE,
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) { 
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
    
  # Checking non-discretionary inputs/outputs
  if ((!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
    warning("This model does not take into account the non-discretionary feature for inputs/outputs.")
  }
    
  # Checking undesirable inputs/outputs
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
      
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
    
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
    
  input <- datadea$input
  output <- datadea$output
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  nc_inputs <- datadea$nc_inputs
  nc_outputs <- datadea$nc_outputs
  nnci <- length(nc_inputs)
  nnco <- length(nc_outputs)
  
  # Checking weights
  if(is.null(weight_slack_i)){
    weight_slack_i <- 1
  }
  if(is.null(weight_slack_o)){
    weight_slack_o <- 1
  }
  
  if (is.matrix(weight_slack_i)) {
    if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
      stop("Invalid weight input matrix (number of inputs x number of evaluated DMUs).")
    }
  } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
    weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
  } else {
    stop("Invalid weight input vector (number of inputs).")
  }
  if ((!is.null(orientation)) && (orientation == "oo")) {
    weight_slack_i <- matrix(0, nrow = ni, ncol = nde)
  }
  rownames(weight_slack_i) <- inputnames
  colnames(weight_slack_i) <- dmunames[dmu_eval]
  
  if (is.matrix(weight_slack_o)) {
    if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
      stop("Invalid weight output matrix (number of outputs x number of evaluated DMUs).")
    }
  } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
    weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
  } else {
    stop("Invalid weight output vector (number of outputs).")
  }
  if ((!is.null(orientation)) && (orientation == "io")) {
    weight_slack_o <- matrix(0, nrow = no, ncol = nde)
  }
  rownames(weight_slack_o) <- outputnames
  colnames(weight_slack_o) <- dmunames[dmu_eval]
  
  # Checking method
  method <- tolower(method)
  method <- match.arg(method)
  
  target_input <- NULL
  target_output <- NULL
    
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  if (method == "milp") { ########################################## MILP method
    
    # Compute extreme efficient DMUs
    if (is.null(extreff)) {
      extreff <- extreme_efficient(datadea = datadea,
                                   dmu_ref = dmu_ref,
                                   rts = rts,
                                   tol = tol)
    } else {
      if (!all(extreff %in% dmu_ref)) {
        stop("Invalid set of extreme efficient DMUs (extreff). It must be a subset of dmu_ref.")
      }
    }
    
    result_add <- model_additive(datadea = datadea,
                                 dmu_eval = extreff,
                                 dmu_ref = dmu_ref,
                                 rts = rts)
    for (jj in 1:length(extreff)) {
      i <- which(dmu_eval == extreff[jj])
      if (length(i) == 1) DMU[[i]] <- result_add$DMU[[jj]]
    }
    
    xineffDMUs <- dmu_eval[!dmu_eval %in% extreff]
    nxineffDMUs <- length(xineffDMUs)
    
    if (nxineffDMUs > 0) {
      
      ndr <- length(extreff)
      inputref <- matrix(input[, extreff], nrow = ni) 
      outputref <- matrix(output[, extreff], nrow = no)
      
      if (rts == "crs") {
        f.con.rs <- NULL
        f.dir.rs <- NULL
        f.rhs.rs <- NULL
        modify <- FALSE
      } else {
        modify <- TRUE
        f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2*ni + 2*no + 3*ndr))
        f.rhs.rs <- 1
        if (rts == "vrs") {
          f.dir.rs <- "="
          M_lambda <- 1
        } else if (rts == "nirs") {
          f.dir.rs <- "<="
          M_lambda <- 1
        } else {
          f.dir.rs <- ">="
        }
      }
      
      # M_d value
      if (is.null(M_d)) {
        if (modify) {
          M_d <- max(inputref)
        } else {
          M_d <- 1e3 * max(inputref)
        }
      }
      
      # Constraints matrix
      f.con.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = ni + 2*no + 3*ndr))
      f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no),
                       matrix(0, nrow = no, ncol = ni + no + 3*ndr))
      f.con.nc <- matrix(0, nrow = (nnci + nnco), ncol = 4*ndr + 2*ni + 2*no)
      f.con.nc[, ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
      f.con.3 <- cbind(matrix(0, nrow = ndr, ncol = ndr + ni + no), -t(inputref),
                       t(outputref), diag(ndr), matrix(0, nrow = ndr, ncol = 2*ndr))
      f.con.nob <- cbind(matrix(0, nrow = ndr, ncol = 2*ndr + 2*ni + 2*no),
                         diag(ndr), diag(ndr))
      f.con.nu <- cbind(matrix(0, nrow = ni, ncol = ndr + ni + no), diag(ni),
                        matrix(0, nrow = ni, ncol = 3*ndr + no))
      f.con.mu <- cbind(matrix(0, nrow = no, ncol = ndr + 2*ni + no), diag(no),
                        matrix(0, nrow = no, ncol = 3*ndr))
      f.con.db <- cbind(matrix(0, nrow = ndr, ncol = ndr + 2*ni + 2*no),
                        -diag(ndr), diag(M_d, ndr), matrix(0, nrow = ndr, ncol = ndr))
      f.con.lamb <- cbind(-diag(ndr), matrix(0, nrow = ndr, ncol = 2*ndr + 2*ni + 2*no),
                          diag(M_lambda, ndr))
      f.con <- rbind(f.con.1, f.con.2, f.con.nc, f.con.3, f.con.nob, f.con.nu,
                     f.con.mu, f.con.db, f.con.lamb, f.con.rs)
      
      # Directions vector
      f.dir <- c(rep("=", ni + no + nnci + nnco + 2*ndr),
                 rep(">=", ni + no + 2*ndr), f.dir.rs)
      
      # Binary variables vector
      binary.vec <- (2*ndr + 2*ni + 2*no + 1) : (4*ndr + 2*ni + 2*no)
      
      if (modify) {
        f.con.db <- cbind(matrix(0, nrow = ndr, ncol = ndr + 2*ni + 2*no),
                          -diag(ndr), diag(M_d, ndr), matrix(0, nrow = ndr, ncol = ndr))
        f.con.aux <- cbind(matrix(0, nrow = 1, ncol = ndr + ni + no),
                           matrix(1, nrow = 1, ncol = ni + no),
                           matrix(0, nrow = 1, ncol = 3*ndr))
        f.con <- rbind(f.con.1, f.con.2, f.con.nc, f.con.3, f.con.nob, f.con.nu,
                       f.con.mu, f.con.db, f.con.lamb, f.con.aux, f.con.rs)
        f.con <- cbind(f.con, matrix(0, nrow = nrow(f.con), ncol = 2))
        f.con[(ni + no + nnci + nnco + 1):(ni + no + nnci + nnco + ndr), ncol(f.con) - 1] <- 1
        f.con[(ni + no + nnci + nnco + 1):(ni + no + nnci + nnco + ndr), ncol(f.con)] <- -1
        f.dir <- c(rep("=", ni + no + nnci + nnco + 2*ndr),
                   rep(">=", ni + no + 2*ndr), "=", f.dir.rs)
      }
      
      for (j in 1:nxineffDMUs) {
        
        i <- which(dmu_eval == xineffDMUs[j])
        ii <- dmu_eval[i]
        
        # Objective function coefficients
        f.obj <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i], rep(0, 3*ndr + ni + no))
        
        # Right hand side vector
        f.rhs <- c(input[, ii], output[, ii], rep(0, nnci + nnco + ndr), rep(1, ndr),
                   weight_slack_i[, i], weight_slack_o[, i], rep(0, 2*ndr), f.rhs.rs)
        #rep(1, ni + no), rep(0, 2*ndr), f.rhs.rs)
        
        if (modify) {
          f.obj <- c(f.obj, 0, 0)
          f.rhs <- c(input[, ii], output[, ii], rep(0, nnci + nnco + ndr), rep(1, ndr),
                     tol*weight_slack_i[, i], tol*weight_slack_o[, i], rep(0, 2*ndr), 1, f.rhs.rs)
        }
        
        if (returnlp) {
          
          lambda <- rep(0, ndr)
          names(lambda) <- dmunames[extreff]
          slack_input <- rep(0, ni)
          names(slack_input) <- inputnames
          slack_output <- rep(0, no)
          names(slack_output) <- outputnames
          var_rest <- NULL
          if (modify) {
            var_rest <- c(0, 0)
          }
          var <- list(lambda = lambda, slack_input = slack_input, slack_output = slack_output,
                      nu = slack_input, mu = slack_output, aux_d = lambda, aux_b = lambda, aux_nob = lambda,
                      var_rest = var_rest)
          DMU[[i]] <- list(direction = "min", objective.in = f.obj, const.mat = f.con,
                           const.dir = f.dir, const.rhs = f.rhs, var = var, binary.vec = binary.vec)
          
        } else {
          
          # if (Rsym){
          #   types <- rep("C", 4*ndr + 2*ni + 2*no)
          #   types[binary.vec] <- "B"
          #   f.dir[which(f.dir == "=")] <- "=="
          #   res <- Rsymphony_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs,
          #                             max = FALSE, types = types)
          # } else {
          res <- lp("min", f.obj, f.con, f.dir, f.rhs, binary.vec = binary.vec)
          # }
          
          if (res$status == 0) {
            
            objval <- res$objval
            
            lambda <- res$solution[1 : ndr]
            names(lambda) <- dmunames[extreff]
            
            slack_input <- res$solution[(ndr + 1) : (ndr + ni)]
            names(slack_input) <- inputnames
            slack_output <- res$solution[(ndr + ni + 1) : (ndr + ni + no)]
            names(slack_output) <- outputnames
            nu <- res$solution[(ndr + ni + no + 1) : (ndr + 2*ni + no)]
            names(nu) <- inputnames
            mu <- res$solution[(ndr + 2*ni + no + 1) : (ndr + 2*ni + 2*no)]
            names(mu) <- outputnames
            aux_d <- res$solution[(ndr + 2*ni + 2*no + 1) : (2*ndr + 2*ni + 2*no)]
            names(aux_d) <- dmunames[extreff]
            aux_b <- res$solution[(2*ndr + 2*ni + 2*no + 1) : (3*ndr + 2*ni + 2*no)]
            names(aux_b) <- dmunames[extreff]
            aux_nob <- res$solution[(3*ndr + 2*ni + 2*no + 1) : (4*ndr + 2*ni + 2*no)]
            names(aux_nob) <- dmunames[extreff]
            var_rest <- NULL
            if (modify) {
              var_rest <- res$solution[(4*ndr + 2*ni + 2*no + 1) : (4*ndr + 2*ni + 2*no + 2)]
            }
            
            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              target_output <- as.vector(outputref %*% lambda)
              names(target_input) <- inputnames
              names(target_output) <- outputnames
            }
            
          } else {
            
            objval <- NA
            lambda <- NA
            slack_input <- NA
            slack_output <- NA
            nu <- NA
            mu <- NA
            aux_d <- NA
            aux_b <- NA
            aux_nob <- NA
            var_rest <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
            }
          }
          
          DMU[[i]] <- list(objval = objval,
                           lambda = lambda,
                           slack_input = slack_input, slack_output = slack_output,
                           nu = nu, mu = mu, aux_d = aux_d, aux_b = aux_b, aux_nob = aux_nob, var_rest = var_rest,
                           target_input = target_input, target_output = target_output)
          
        }
      }
    }
    
  } else { ############################################## Maximal friends method 
    
    # Compute maximal friends
    if (is.null(maxfr)) {
      if (!silent) {
        print("Computing maximal friends...")
      }
      maxfr <- maximal_friends(datadea = datadea,
                               dmu_ref = dmu_ref,
                               rts = rts,
                               tol = tol,
                               silent = silent)
    } else {
      if (!all(unique(unlist(maxfr)) %in% dmu_ref)) {
        stop("Invalid set of maximal friends (maxfr). They must be subsets of dmu_ref.")
      }
    }
    
    effDMUs <- sort(unique(unlist(maxfr)))
    result_add <- model_additive(datadea = datadea,
                                 dmu_eval = effDMUs,
                                 dmu_ref = dmu_ref,
                                 rts = rts)
    for (jj in 1:length(effDMUs)) {
      i <- which(dmu_eval == effDMUs[jj])
      if (length(i) > 0) DMU[[i]] <- result_add$DMU[[jj]]
    }
    
    ineffDMUs <- dmu_eval[!dmu_eval %in% effDMUs]
    nineffDMUs <- length(ineffDMUs)
    
    if (nineffDMUs > 0) {
      
      nmf <- length(maxfr)
      all_kaizen <- vector(mode = "list", length = nineffDMUs)
      names(all_kaizen) <- dmunames[ineffDMUs]
      dmu_ref_orig <- dmu_ref
      
      for (j in 1:nineffDMUs) {
        
        i <- which(dmu_eval == ineffDMUs[j])
        ii <- dmu_eval[i]
        all_kaizen[[j]] <- vector(mode = "list", length = nmf)
        names(all_kaizen[[j]]) <- names(maxfr)
        
        for (k in 1:nmf) {
          
          dmu_ref <- maxfr[[k]]
          ndr <- length(dmu_ref)
          inputref <- matrix(input[, dmu_ref], nrow = ni) 
          outputref <- matrix(output[, dmu_ref], nrow = no)
          
          if (rts == "crs") {
            f.con.rs <- NULL
            f.dir.rs <- NULL
            f.rhs.rs <- NULL
          } else {
            f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
            f.rhs.rs <- 1
            if (rts == "vrs") {
              f.dir.rs <- "="
            } else if (rts == "nirs") {
              f.dir.rs <- "<="
            } else {
              f.dir.rs <- ">="
            }
          }
          
          # Directions vector
          f.dir <- c(rep("=", ni + no + nnci + nnco), f.dir.rs)
          
          # Right hand side vector
          f.rhs <- c(input[, ii], output[, ii], rep(0, nnci + nnco), f.rhs.rs)
          
          # Objective function coefficients
          f.obj <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i])
          
          # Constraints matrix
          f.con.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = no))
          f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no))
          f.con.nc <- matrix(0, nrow = nnci + nnco, ncol = ndr + ni + no)
          f.con.nc[, ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
          f.con <- rbind(f.con.1, f.con.2, f.con.nc, f.con.rs)
          
          res <- lp("min", f.obj, f.con, f.dir, f.rhs)
          
          if (res$status == 0) {
            
            objval <- res$objval
            
            lambda <- res$solution[1 : ndr]
            names(lambda) <- dmunames[dmu_ref]
            
            slack_input <- res$solution[(ndr + 1) : (ndr + ni)]
            names(slack_input) <- inputnames
            slack_output <- res$solution[(ndr + ni + 1) : (ndr + ni + no)]
            names(slack_output) <- outputnames
            
            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              target_output <- as.vector(outputref %*% lambda)
              names(target_input) <- inputnames
              names(target_output) <- outputnames
            }
            
          } else {
            
            objval <- 1e12
            lambda <- NA
            slack_input <- NA
            slack_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
            }
            
          }
          
          all_kaizen[[j]][[k]] <- list(objval = objval,
                                       lambda = lambda,
                                       slack_input = slack_input, slack_output = slack_output,
                                       target_input = target_input, target_output = target_output)
          
        }
      }
      
      dmu_ref <- dmu_ref_orig
      
      for (j in 1:nineffDMUs) {
        kbest <- 1
        if (nmf > 1) {
          for (k in 2:nmf) {
            if (all_kaizen[[j]][[k]]$objval < all_kaizen[[j]][[kbest]]$objval) {
              kbest <- k
            }
          }
        }
        i <- which(dmu_eval == ineffDMUs[j])
        DMU[[i]] <- all_kaizen[[j]][[kbest]]
      }
    }
  }
  
  if (check_target) {
    target_input_list <- lapply(DMU, function(x) x$target_input)
    target_output_list <- lapply(DMU, function(x) x$target_output)
    for (i in 1:nde) {
      ii <- dmu_eval[i]
      datadea2 <- datadea
      datadea2$input[, ii] <- target_input_list[[i]]
      datadea2$output[, ii] <- target_output_list[[i]]
      objval <- model_additive(datadea = datadea2, dmu_eval = i,
                               dmu_ref = dmu_ref, rts = rts)$DMU[[1]]$objval
      if (objval > tol * (max(target_input_list[[i]], target_output_list[[i]]))) {
        warning(paste("DMU", ii, "named", dmunames[ii], "has inefficient target. Additive objval: ", objval))
      }
    }
  }
  
  deaOutput <- list(modelname = "addmin",
                    rts = rts,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    weight_slack_i = weight_slack_i,
                    weight_slack_o = weight_slack_o,
                    orientation = NA)
  
  return(structure(deaOutput, class = "dea"))
  
}