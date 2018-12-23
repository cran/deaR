#' @title Additive super-efficiency DEA model.
#'   
#' @description Solve the additive super-efficiency model proposed by Du, Liang and Zhu (2010). It is an extension of the SBM super-efficiency to the additive DEA model.
#' 
#' @usage model_addsupereff(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             weight_slack_i = NULL,
#'             weight_slack_o = NULL,
#'             rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'             L = 1,
#'             U = 1,
#'             compute_target = TRUE,
#'             returnlp = FALSE,
#'             ...)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x \code{ne} (where \code{ne} is the lenght of \code{dmu_eval})
#'                   with the weights of the input superslacks (\code{t_input}).
#'                   If \code{weight_slack_i} is the matrix of the inverses of inputs (of DMUS in \code{dmu_eval}), the model is unit invariant.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x \code{ne} (where \code{ne} is the lenght of \code{dmu_eval})
#'                   with the weights of the output superslacks (\code{t_output}).
#'                   If \code{weight_slack_o} is the matrix of the inverses of outputs (of DMUS in \code{dmu_eval}), the model is unit invariant.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#'            "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets, projections and slacks. 
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective function and constraints).
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
#' Du, J.; Liang, L.; Zhu, J. (2010). "A Slacks-based Measure of Super-efficiency in Data Envelopment Analysis. A Comment", European Journal of Operational Research, 204, 694-697. \url{https://doi.org/10.1016/j.ejor.2009.12.007}
#' 
#' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking. Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York. Doi: 10.1007/978-3-319-06647-9.
#' 
#' @examples 
#' # Replication of results in Du, Liang and Zhu (2010, Table 6, p.696)
#' data("Power_plants")
#' Power_plants <- read_data(Power_plants, 
#'                           ni = 4, 
#'                           no = 2)
#' result <- model_addsupereff(Power_plants, 
#'                             rts = "crs")  
#' efficiencies(result)
#' 
#' @seealso \code{\link{model_additive}}, \code{\link{model_supereff}},  \code{\link{model_sbmsupereff}}
#' 
#' @import lpSolve
#' 
#' @export

model_addsupereff <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           weight_slack_i = NULL,
           weight_slack_o = NULL,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           compute_target = TRUE,
           returnlp = FALSE,
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
  #Checking non-controllable or non-discretionary inputs/outputs
  #if ((!is.null(datadea$nc_inputs)) || (!is.null(datadea$nc_outputs))
  #    || (!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
  #  warning("This model does not take into account non-controllable or non-discretionary feature for inputs/outputs.")
  #  datadea$nc_inputs <- NULL
  #  datadea$nc_outputs <- NULL
  #  datadea$nd_inputs <- NULL
  #  datadea$nd_outputs <- NULL
  #}
      
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
    
  input <- datadea$input
  output <- datadea$output
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  # Zeros in input and output data. Case 2 (Tone 2001)
  nzimin <- apply(input, MARGIN = 1, function(x) min(x[x > 0])) / 100
  nzomin <- apply(output, MARGIN = 1, function(x) min(x[x > 0])) / 100
  for (ii in dmu_eval) {
    input[which(input[, ii] == 0), ii] <- nzimin[which(input[, ii] == 0)]
    output[which(output[, ii] == 0), ii] <- nzomin[which(output[, ii] == 0)]
  }
  
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  nc_inputs <- datadea$nc_inputs
  nc_outputs <- datadea$nc_outputs
  nnci <- length(nc_inputs)
  nnco <- length(nc_outputs)
  
  # Checking weights
  if (is.null(weight_slack_i)) {
    weight_slack_i <- matrix(1 / input[, dmu_eval], nrow = ni) / (ni + no - nnci - nnco)
  } else {
    if (is.matrix(weight_slack_i)) {
      if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
        stop("Invalid matrix of weights of the input slacks (number of inputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
      weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
    } else {
      stop("Invalid vector of weights of the input slacks.")
    }
  }
  weight_slack_i[nc_inputs, ] <- 0
  rownames(weight_slack_i) <- inputnames
  colnames(weight_slack_i) <- dmunames[dmu_eval]
  
  if (is.null(weight_slack_o)) {
    weight_slack_o <- matrix(1 / output[, dmu_eval], nrow = no) / (ni + no - nnci - nnco)
  } else {
    if (is.matrix(weight_slack_o)) {
      if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
        stop("Invalid matrix of weights of the output slacks (number of outputs x number of evaluated DMUs).")
      }
    } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
      weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
    } else {
      stop("Invalid vector of weights of the output slacks.")
    }
  }
  weight_slack_o[nc_outputs, ] <- 0
  rownames(weight_slack_o) <- outputnames
  colnames(weight_slack_o) <- dmunames[dmu_eval]
  
  target_input <- NULL
  target_output <- NULL
  project_input <- NULL
  project_output <- NULL
  slack_input <- NULL
  slack_output <- NULL
    
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
    
  ###########################
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
    if (rts == "vrs") {
      f.dir.rs <- "="
      f.rhs.rs <- 1
    } else if (rts == "nirs") {
      f.dir.rs <- "<="
      f.rhs.rs <- 1
    } else if (rts == "ndrs") {
      f.dir.rs <- ">="
      f.rhs.rs <- 1
    } else {
      f.con.rs <- rbind(f.con.rs, f.con.rs)
      f.dir.rs <- c(">=", "<=")
      f.rhs.rs <- c(L, U)
    }
  }
  
  # Matriz técnica
  f.con.1 <- cbind(inputref, -diag(ni), matrix(0, nrow = ni, ncol = no))
  f.con.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), diag(no))
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    w0i <- which(weight_slack_i[, i] == 0)
    nw0i <- length(w0i)
    w0o <- which(weight_slack_o[, i] == 0)
    nw0o <- length(w0o)
    
    # Vector de coeficientes de la función objetivo
    f.obj <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i])
    
    # Matriz técnica
    f.con.se <- rep(0, ndr)
    f.con.se[dmu_ref == ii] <- 1
    f.con.se <- c(f.con.se, rep(0, ni + no))
    f.con.w0 <- matrix(0, nrow = (nw0i + nw0o), ncol = (ndr + ni + no))
    f.con.w0[, ndr + c(w0i, ni + w0o)] <- diag(nw0i + nw0o)
    f.con <- rbind(f.con.1, f.con.2, f.con.se, f.con.w0, f.con.rs)
    
    # Vector de dirección de restricciones
    f.dir <- c(rep("<=", ni), rep(">=", no), rep("=", 1 + nw0i + nw0o), f.dir.rs)
    f.dir[nc_inputs] <- "="
    f.dir[ni + nc_outputs] <- "="
    
    # Vector de términos independientes
    f.rhs <- c(input[, ii], output[, ii], rep(0, 1 + nw0i + nw0o), f.rhs.rs)
    
    if (returnlp) {
      
      lambda <- rep(0, ndr)
      names(lambda) <- dmunames[dmu_ref]
      t_input <- rep(0, ni)
      names(t_input) <- inputnames
      t_output <- rep(0, no)
      names(t_output) <- outputnames
      var <- list(lambda = lambda, t_input = t_input, t_output = t_output)
      DMU[[i]] <- list(direction = "min", objective.in = f.obj, const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs,
                       var = var)
      
    } else {
    
      res <- lp("min", f.obj, f.con, f.dir, f.rhs)
      
      if (res$status == 0) {
        
        objval <- res$objval
        
        lambda <- res$solution[1 : ndr]
        names(lambda) <- dmunames[dmu_ref]
        
        t_input <- res$solution[(ndr + 1) : (ndr + ni)] # input superslack
        names(t_input) <- inputnames
        t_output <- res$solution[(ndr + ni + 1) : (ndr + ni + no)] # output superslack
        names(t_output) <- outputnames
        
        delta_num <- 1 + sum(t_input / input[, ii]) / (ni - nnci)
        delta_den <- 1 - sum(t_output / output[, ii]) / (no - nnco)
        delta <- delta_num / delta_den
        
        if (compute_target) {
          target_input <- as.vector(inputref %*% lambda)
          names(target_input) <- inputnames
          target_output <- as.vector(outputref %*% lambda)
          names(target_output) <- outputnames
          
          project_input <- input[, ii] + t_input
          names(project_input) <- inputnames
          project_output <- output[, ii] - t_output
          names(project_output) <- outputnames
          
          slack_input <- project_input - target_input
          names(slack_input) <- inputnames
          slack_output <- target_output - project_output
          names(slack_output) <- outputnames
        }
        
      } else {
        
        delta <- NA
        objval <- NA
        lambda <- NA
        t_input <- NA
        t_output <- NA
        if (compute_target) {
          target_input <- NA
          target_output <- NA
          project_input <- NA
          project_output <- NA
          slack_input <- NA
          slack_output <- NA
        }
        
      }

      DMU[[i]] <- list(delta = delta,
                       objval = objval,
                       lambda = lambda,
                       t_input = t_input, t_output = t_output,
                       slack_input = slack_input, slack_output = slack_output,
                       project_input = project_input, project_output = project_output,
                       target_input = target_input, target_output = target_output)
      
    }
    
  }
  
  deaOutput <- list(modelname = "addsupereff",
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    weight_slack_i = weight_slack_i,
                    weight_slack_o = weight_slack_o,
                    orientation = "NA")
  
  return(structure(deaOutput, class = "dea"))
  
}