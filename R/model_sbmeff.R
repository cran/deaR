#' @title Slack based measure (SBM) of efficiency model.
#'   
#' @description Calculate the SBM model proposed by Tone (2001).
#' 
#' @usage model_sbmeff(datadea,
#'              dmu_eval = NULL,
#'              dmu_ref = NULL,
#'              weight_input = 1,
#'              weight_output = 1,
#'              orientation = c("no", "io", "oo"),
#'              rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'              L = 1,
#'              U = 1,
#'              kaizen = FALSE,
#'              maxfr = NULL,
#'              tol = 1e-6,
#'              silent = FALSE,
#'              compute_target = TRUE,
#'              returnlp = FALSE,
#'              ...)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param weight_input A value, vector of length \code{m}, or matrix \code{m} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval}) with weights to
#' inputs corresponding to the relative importance of items.
#' @param weight_output A value, vector of length \code{m}, or matrix \code{m} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval}) with weights to
#' outputs corresponding to the relative importance of items.
#' @param orientation A string, equal to "no" (non-oriented), "io" (input-oriented)
#' or "oo" (output-oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param kaizen Logical. If \code{TRUE}, the kaizen version of SBM (Tone 2010),
#' also known as SBM-Max, is computed.
#' @param maxfr A list with the maximal friends sets, as it is returned by function
#' \code{maximal_friends}. If \code{NULL} (default) this list is computed internally.
#' @param tol Numeric, a tolerance margin for checking efficiency (only for the kaizen version).
#' @param silent Logical. If \code{FALSE} (default) it prints all the messages from
#' function \code{maximal_friends}.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets. 
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems
#' (objective function and constraints). If \code{kaizen} is \code{TRUE} it is ignored.
#' @param ... Other options (currently not implemented)
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
#' Tone, K. (2001). "A Slacks-Based Measure of Efficiency in Data Envelopment Analysis",
#' European Journal of Operational Research, 130, 498-509. \doi{10.1016/S0377-2217(99)00407-5}
#'
#' Tone, K. (2010). "Variations on the theme of slacks-based measure of efficiency in DEA",
#' European Journal of Operational Research, 200, 901-907. \doi{10.1016/j.ejor.2009.01.027}
#'
#' Cooper, W.W.; Seiford, L.M.; Tone, K. (2007). Data Envelopment Analysis. A Comprehensive
#' Text with Models, Applications, References and DEA-Solver Software. 2nd Edition. Springer,
#' New York. \doi{10.1007/978-0-387-45283-8}
#' 
#' @examples 
#' # Replication of results in Tone (2001, p.505)
#' data("Tone2001")
#' data_example <- read_data(Tone2001, 
#'                           ni = 2, 
#'                           no = 2)
#' result_SBM <- model_sbmeff(data_example, 
#'                            orientation = "no", 
#'                            rts = "crs")
#' result_CCR <- model_basic(data_example, 
#'                           orientation = "io", 
#'                           rts = "crs")
#' efficiencies(result_SBM)
#' efficiencies(result_CCR)
#' slacks(result_SBM)
#' slacks(result_CCR)
#'  
#' # Example. Replication of results in Tone (2003), pp 10-11 case 1:1.
#' data("Tone2003")
#' data_example <- read_data(Tone2003,
#'                           ni = 1,
#'                           no = 2,
#'                           ud_outputs = 2)
#' result <- model_sbmeff(data_example,
#'                        rts = "vrs")
#' efficiencies(result)
#' targets(result)
#'  
#'  
#' @seealso \code{\link{model_nonradial}}, \code{\link{model_deaps}},
#' \code{\link{model_profit}}, \code{\link{model_sbmsupereff}}
#' 
#' @import lpSolve
#' 
#' @export

model_sbmeff <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           weight_input = 1,
           weight_output = 1,
           orientation = c("no", "io", "oo"),
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           kaizen = FALSE,
           maxfr = NULL,
           tol = 1e-6,
           silent = FALSE,
           compute_target = TRUE,
           returnlp = FALSE, ...) {
    
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
    
  input <- datadea$input
  output <- datadea$output
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  # Zeros in output data. Case 2 (Tone 2001)
  nzomin <- apply(output, MARGIN = 1, function(x) min(x[x > 0])) / 100
  for (ii in dmu_eval) {
    output[which(output[, ii] == 0), ii] <- nzomin[which(output[, ii] == 0)]
  }
  
  nc_inputs <- datadea$nc_inputs
  nc_outputs <- datadea$nc_outputs
  nnci <- length(nc_inputs)
  nnco <- length(nc_outputs)
  ud_inputs <- datadea$ud_inputs
  if ((is.null(ud_inputs) == FALSE) && orientation != "oo"){
    warning("Undesirable inputs with not output-oriented model could generate
            negative efficiencies. The lower the efficiency, the more inefficient the DMU.")
  }
  ud_outputs <- datadea$ud_outputs
  aux_udi <- rep(1, ni)
  aux_udi[ud_inputs] <- -1
  aux_udo <- rep(1, no)
  aux_udo[ud_outputs] <- -1
  
  aux_i <- 1
  aux_o <- 1
  if (orientation == "io") {
    aux_o <- 0
  } else if (orientation == "oo") {
    aux_i <- 0
  }
  
  # Checking weights
  if (is.matrix(weight_input)) {
    if ((nrow(weight_input) != ni) || (ncol(weight_input) != nde)) {
      stop("Invalid input weights matrix (number of inputs x number of evaluated DMUs).")
    }
  } else if ((length(weight_input) == 1) || (length(weight_input) == ni)) {
    weight_input <- matrix(weight_input, nrow = ni, ncol = nde)
  } else {
    stop("Invalid input weights vector (number of inputs).")
  }
  weight_input[nc_inputs, ] <- 0
  sumwi <- colSums(weight_input)
  if (any(sumwi == 0) && aux_i == 1) {
    stop("A sum of input weights is 0.")
  }
  rownames(weight_input) <- inputnames
  colnames(weight_input) <- dmunames[dmu_eval]
  
  if (is.matrix(weight_output)) {
    if ((nrow(weight_output) != no) || (ncol(weight_output) != nde)) {
      stop("Invalid output weights matrix (number of outputs x number of evaluated DMUs).")
    }
  } else if ((length(weight_output) == 1) || (length(weight_output) == no)) {
    weight_output <- matrix(weight_output, nrow = no, ncol = nde)
  } else {
    stop("Invalid output weights vector (number of outputs).")
  }
  weight_output[nc_outputs, ] <- 0
  sumwo <- colSums(weight_output)
  if (any(sumwo == 0) && aux_o == 1) {
    stop("A sum of output weights is 0.")
  }
  rownames(weight_output) <- outputnames
  colnames(weight_output) <- dmunames[dmu_eval]
  
  target_input <- NULL
  target_output <- NULL
  all_kaizen <- NULL
    
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  if (!kaizen) { #################### Compute SBM-Min efficiency (Tone 2001) ####################
    
    ndr <- length(dmu_ref)
    inputref <- matrix(input[, dmu_ref], nrow = ni) 
    outputref <- matrix(output[, dmu_ref], nrow = no)
    
    if (rts == "crs") {
      f.con.rs <- NULL
      f.dir.rs <- NULL
      f.rhs.rs <- NULL
    } else {
      f.con.rs <- cbind(-1, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
      f.rhs.rs <- 0
      if (rts == "vrs") {
        f.dir.rs <- "="
      } else if (rts == "nirs") {
        f.dir.rs <- "<="
      } else if (rts == "ndrs") {
        f.dir.rs <- ">="
      } else {
        f.con.rs <- rbind(cbind(-L, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no)),
                          cbind(-U, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no)))
        f.dir.rs <- c(">=", "<=")
        f.rhs.rs <- c(0, 0)
      }
    }
    
    # Constraints matrix
    f.con.nc <- matrix(0, nrow = (nnci + nnco), ncol = (1 + ndr + ni + no))
    f.con.nc[, 1 + ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
    
    # Directions vector
    f.dir <- c(rep("=", 1 + ni + no + nnci + nnco), f.dir.rs)
    
    # Right hand side vector
    f.rhs <- c(1, rep(0, ni + no + nnci + nnco), f.rhs.rs)
    
    for (i in 1:nde) {
      
      ii <- dmu_eval[i]
      
      # Zeros in input data
      zero_inputs <- which(input[, ii] == 0)
      nzinput <- input[, ii]
      nzinput[zero_inputs] <- 1 # zero inputs become 1
      weight_input[zero_inputs, i] <- 0 # weights corresponding to zero inputs become 0
      nzi_sumwi <- 1
      if (aux_i == 1) {
        nzi_sumwi <- sum(weight_input[, i])
      }
      if (nzi_sumwi == 0) {
          stop("A sum of nonzero-input weights is 0.")
      }
      
      # Objective function coefficients
      f.obj <- c(1, rep(0, ndr), -aux_i * weight_input[, i] / (nzi_sumwi * nzinput), rep(0, no))
      
      # Constraints matrix
      f.con.0 <- c(1, rep(0, ndr + ni), aux_o * weight_output[, i] / (sumwo[i] * output[, ii]))
      f.con.1 <- cbind(-input[, ii], inputref, diag(aux_udi), matrix(0, nrow = ni, ncol = no))
      f.con.2 <- cbind(-output[, ii], outputref, matrix(0, nrow = no, ncol = ni), -diag(aux_udo))
      f.con <- rbind(f.con.0, f.con.1, f.con.2, f.con.nc, f.con.rs)
      
      if (returnlp) {
        
        t <- 0
        names(t) <- "t"
        tlambda <- rep(0, ndr)
        names(tlambda) <- dmunames[dmu_ref]
        tslack_input <- rep(0, ni)
        names(tslack_input) <- inputnames
        tslack_output <- rep(0, no)
        names(tslack_output) <- outputnames
        var <- list(t = t, tlambda = tlambda, tslack_input = tslack_input,
                    tslack_output = tslack_output)
        DMU[[i]] <- list(direction = "min", objective.in = f.obj, const.mat = f.con,
                         const.dir = f.dir, const.rhs = f.rhs, var = var)
        
      } else {
        
        res <- lp("min", f.obj, f.con, f.dir, f.rhs)
        
        if (res$status == 0) {
          
          efficiency <- res$objval
          res <- res$solution
          
          t <- res[1]
          lambda <- res[2 : (ndr + 1)] / t
          names(lambda) <- dmunames[dmu_ref]
          
          slack_input <- res[(ndr + 2) : (ndr + ni + 1)] / t
          names(slack_input) <- inputnames
          slack_output <- res[(ndr + ni + 2) : (ndr + ni + no + 1)] / t
          names(slack_output) <- outputnames
          
          if (compute_target) {
            target_input <- as.vector(inputref %*% lambda)
            names(target_input) <- inputnames
            target_output <- as.vector(outputref %*% lambda)
            names(target_output) <- outputnames
          }
          
        } else {

          efficiency <- NA
          lambda <- NA
          slack_input <- NA
          slack_output <- NA
          if (compute_target) {
            target_input <- NA
            target_output <- NA
          }
          
        }
        
        DMU[[i]] <- list(efficiency = efficiency,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output)
        
      }
      
    }
    
  } else { #################### Compute Kaizen efficiency (Tone 2010) ####################
    
    if (orientation != "no") {
    #  orientation <- "no"
    #  warning("Orientation changed to non-oriented.")
      warning("Orientation is not non-oriented.")
    }
    
    if (rts == "grs") {
      stop("Kaizen is not available for generalized returns to scale.")
    }
    
    # Find efficient/inefficient DMUs in dmu_eval
    result_sbm <- model_sbmeff(datadea = datadea,
                               dmu_eval = dmu_eval,
                               dmu_ref = dmu_ref,
                               rts = rts)
    eff_sbm <- efficiencies(result_sbm)
    effDMUs <- dmu_eval[eff_sbm >= (1 - tol)]
    ineffDMUs <- dmu_eval[!dmu_eval %in% effDMUs]
    nineffDMUs <- length(ineffDMUs)
    
    for (jj in effDMUs) {
      i <- which(dmu_eval == jj)
      DMU[[i]] <- result_sbm$DMU[[i]]
    }
    
    if (nineffDMUs > 0) {
      
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
      
      nmf <- length(maxfr)
      all_kaizen <- vector(mode = "list", length = nineffDMUs)
      names(all_kaizen) <- dmunames[ineffDMUs]
      dmu_ref_orig <- dmu_ref
      
      for (j in 1:nineffDMUs) {
        
        i <- which(dmu_eval == ineffDMUs[j])
        ii <- dmu_eval[i]
        all_kaizen[[j]] <- vector(mode = "list", length = nmf)
        names(all_kaizen[[j]]) <- names(maxfr)
        
        # Zeros in input data
        zero_inputs <- which(input[, ii] == 0)
        nzinput <- input[, ii]
        nzinput[zero_inputs] <- 1 # zero inputs become 1
        weight_input[zero_inputs, i] <- 0 # weights corresponding to zero inputs become 0
        nzi_sumwi <- 1
        if (aux_i == 1) {
          nzi_sumwi <- sum(weight_input[, i])
        }
        if (nzi_sumwi == 0) {
          stop("A sum of nonzero-input weights is 0.")
        }
        
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
            f.con.rs <- cbind(-1, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
            f.rhs.rs <- 0
            if (rts == "vrs") {
              f.dir.rs <- "="
            } else if (rts == "nirs") {
              f.dir.rs <- "<="
            } else {
              f.dir.rs <- ">="
            }
          }
          
          # Constraints matrix
          f.con.nc <- matrix(0, nrow = (nnci + nnco), ncol = (1 + ndr + ni + no))
          f.con.nc[, 1 + ndr + c(nc_inputs, ni + nc_outputs)] <- diag(nnci + nnco)
          
          # Directions vector
          f.dir <- c(rep("=", 1 + ni + no + nnci + nnco), f.dir.rs)
          
          # Right hand side vector
          f.rhs <- c(1, rep(0, ni + no + nnci + nnco), f.rhs.rs)
          
          # Objective function coefficients
          f.obj <- c(1, rep(0, ndr), -aux_i * weight_input[, i] / (nzi_sumwi * nzinput), rep(0, no))
          
          # Constraints matrix
          f.con.0 <- c(1, rep(0, ndr + ni), aux_o * weight_output[, i] / (sumwo[i] * output[, ii]))
          f.con.1 <- cbind(-input[, ii], inputref, diag(aux_udi), matrix(0, nrow = ni, ncol = no))
          f.con.2 <- cbind(-output[, ii], outputref, matrix(0, nrow = no, ncol = ni), -diag(aux_udo))
          f.con <- rbind(f.con.0, f.con.1, f.con.2, f.con.nc, f.con.rs)
          
          res <- lp("max", f.obj, f.con, f.dir, f.rhs)
          
          if (res$status == 0) {
            
            efficiency <- res$objval
            res <- res$solution
            
            t <- res[1]
            lambda <- res[2 : (ndr + 1)] / t
            names(lambda) <- dmunames[dmu_ref]
            
            slack_input <- res[(ndr + 2) : (ndr + ni + 1)] / t
            names(slack_input) <- inputnames
            slack_output <- res[(ndr + ni + 2) : (ndr + ni + no + 1)] / t
            names(slack_output) <- outputnames
            
            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              names(target_input) <- inputnames
              target_output <- as.vector(outputref %*% lambda)
              names(target_output) <- outputnames
            }
            
          } else {
            
            efficiency <- 0
            lambda <- NA
            slack_input <- NA
            slack_output <- NA
            if (compute_target) {
              target_input <- NA
              target_output <- NA
            }
          }
          
          all_kaizen[[j]][[k]] <- list(efficiency = efficiency,
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
            if (all_kaizen[[j]][[k]]$efficiency > all_kaizen[[j]][[kbest]]$efficiency) {
              kbest <- k
            }
          }
        }
        i <- which(dmu_eval == ineffDMUs[j])
        DMU[[i]] <- all_kaizen[[j]][[kbest]]
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
  
  deaOutput <- list(modelname = "sbmeff",
                    orientation = orientation,
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    kaizen = kaizen,
                    maxfr = maxfr,
                    all_kaizen = all_kaizen,
                    weight_input = weight_input,
                    weight_output = weight_output)
  
  return(structure(deaOutput, class = "dea"))
  
  }
