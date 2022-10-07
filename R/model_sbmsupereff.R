#' @title Slack based measure of superefficiency model
#'   
#' @description Slack based measure of superefficiency model (Tone 2002) with \code{n}
#' DMUs, \code{m} inputs and \code{s} outputs.
#' 
#' @usage model_sbmsupereff(datadea,
#'                   dmu_eval = NULL,
#'                   dmu_ref = NULL,
#'                   weight_input = 1,
#'                   weight_output = 1,
#'                   orientation = c("no", "io", "oo"),
#'                   rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                   L = 1,
#'                   U = 1,
#'                   compute_target = TRUE,
#'                   compute_rho = FALSE,
#'                   kaizen = FALSE,
#'                   silent = FALSE,
#'                   returnlp = FALSE)
#' 
#' @param datadea The data, including DMUs, inputs and outputs.
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
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets,
#' superslacks (\code{t_input} and \code{t_output}) and slacks.
#' @param compute_rho Logical. If it is \code{TRUE}, it computes the SBM efficiency
#' score (applying \code{model_sbmeff}) of the DMU (\code{project_input}, \code{project_output}).
#' @param kaizen Logical. If \code{TRUE}, the kaizen version of SBM (Tone 2010),
#' also known as SBM-Max, is computed for the efficiency score of the DMU
#' (\code{project_input}, \code{project_output}).
#' @param silent Logical. If \code{FALSE} (default) it prints all the messages
#' from function \code{maximal_friends}.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems
#' (objective function and constraints).
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
#' Tone, K. (2002). "A slacks-based measure of super-efficiency in data envelopment
#' analysis", European Journal of Operational Research, 143, 32-41.
#' \doi{10.1016/S0377-2217(01)00324-1}  
#'
#' Tone, K. (2010). "Variations on the theme of slacks-based measure of efficiency
#' in DEA", European Journal of Operational Research, 200, 901-907.
#' \doi{10.1016/j.ejor.2009.01.027}
#'   
#' Cooper, W.W.; Seiford, L.M.; Tone, K. (2007). Data Envelopment Analysis.
#' A Comprehensive Text with Models, Applications, References and DEA-Solver Software.
#' 2nd Edition. Springer, New York. \doi{10.1007/978-0-387-45283-8}
#' 
#' @examples 
#' # Replication of results in Tone(2002, p.39)
#' data("Power_plants")
#' data_example <- read_data(Power_plants,
#'                           ni = 4,
#'                           no = 2)
#' result <- model_sbmsupereff(data_example,
#'                             orientation = "io",
#'                             rts = "crs") 
#' efficiencies(result)
#' slacks(result)$slack_input
#' references(result)
#' 
#' @seealso \code{\link{model_sbmeff}}, \code{\link{model_supereff}},
#' \code{\link{model_addsupereff}}
#' 
#' @import lpSolve
#' 
#' @export

model_sbmsupereff <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           weight_input = 1,
           weight_output = 1,
           orientation = c("no", "io", "oo"),
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           compute_target = TRUE,
           compute_rho = FALSE,
           kaizen = FALSE,
           silent = FALSE,
           returnlp = FALSE) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
  # Checking undesirable io and rts
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
  
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
    
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
  # Possible non-feasibilities
  if ((rts != "crs") && (orientation != "no")) {
    warning("For oriented models and non constant returns to scale, feasibility is not
             guaranteed. Proceed with caution, some DMUs results may be missing!")
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
  #ud_inputs <- datadea$ud_inputs
  #ud_outputs <- datadea$ud_outputs
  dir1 <- rep("<=", ni)
  #dir1[ud_inputs] <- ">="
  dir2 <- rep(">=", no)
  #dir2[ud_outputs] <- "<="
  dir3 <- rep(">=", ni)
  #dir3[ud_inputs] <- "<="
  dir4 <- rep("<=", no)
  #dir4[ud_outputs] <- ">="
  
  if (orientation == "oo") {
    obj <- "max"
  } else {
    obj <- "min"
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
  if (any(sumwi == 0)) {
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
  if (any(sumwo == 0)) {
    stop("A sum of output weights is 0.")
  }
  rownames(weight_output) <- outputnames
  colnames(weight_output) <- dmunames[dmu_eval]
  
  target_input <- NULL
  target_output <- NULL
  t_input <- NULL
  t_output <- NULL
  slack_input <- NULL
  slack_output <- NULL
  rho <- NULL
  gamma <- NULL
  csbm <- NULL
    
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
    
  ###########################
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else if (rts == "grs") {
    f.con.rs <- rbind(cbind(-L, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no)),
                      cbind(-U, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no)))
    f.dir.rs <- c(">=", "<=")
    f.rhs.rs <- c(0, 0)
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
  f.con.1 <- cbind(0, inputref, -diag(ni), matrix(0, nrow = ni, ncol = no))
  f.con.2 <- cbind(0, outputref, matrix(0, nrow = no, ncol = ni), -diag(no))
  
  # Directions vector and constraints matrix
  if (orientation == "no") {
    f.dir <- c("=", dir1, dir2, dir3, dir4)
  } else if (orientation == "io") {
    f.con.0 <- c(1, rep(0, ndr + ni + no))
    f.dir <- c("=", dir1, dir2, dir3, rep("=", no))
  } else {
    f.con.0 <- c(1, rep(0, ndr + ni + no))
    f.dir <- c("=", dir1, dir2, rep("=", ni), dir4)
  }
  f.dir[1 + c(nc_inputs, ni + nc_outputs, ni + no + nc_inputs, ni + no + ni + nc_outputs)] <- "="
  f.dir <- c(f.dir, "=", f.dir.rs)
  
  # Right hand side vector
  f.rhs <- c(1, rep(0, ni + no + ni + no + 1), f.rhs.rs)
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    if (orientation == "no") {
      f.obj <- c(0, rep(0, ndr), weight_input[, i] / (sumwi[i] * input[, ii]), rep(0, no))
      f.con.0 <- c(0, rep(0, ndr), rep(0, ni), weight_output[, i] / (sumwo[i] * output[, ii]))
    } else if (orientation == "io") {
      f.obj <- c(0, rep(0, ndr), weight_input[, i] / (sumwi[i] * input[, ii]), rep(0, no))
    } else {
      f.obj <- c(0, rep(0, ndr), rep(0, ni), weight_output[, i] / (sumwo[i] * output[, ii]))
    }
      
    # Constraints matrix
    f.con.3 <- cbind(-input[, ii], matrix(0, nrow = ni, ncol = ndr), diag(ni), matrix(0, nrow = ni, ncol = no))
    f.con.4 <- cbind(-output[, ii], matrix(0, nrow = no, ncol = (ndr + ni)), diag(no))
    f.con.se <- rep(0, ndr)
    f.con.se[dmu_ref == ii] <- 1
    f.con.se <- c(0, f.con.se, rep(0, ni + no))
    f.con <- rbind(f.con.0, f.con.1, f.con.2, f.con.3, f.con.4, f.con.se, f.con.rs)
    
    if (returnlp) {
      
      t <- 0
      names(t) <- "t"
      tlambda <- rep(0, ndr)
      names(tlambda) <- dmunames[dmu_ref]
      tproject_input <- rep(0, ni)
      names(tproject_input) <- inputnames
      tproject_output <- rep(0, no)
      names(tproject_output) <- outputnames
      var <- list(t = t, tlambda = tlambda, tproject_input = tproject_input,
                  tproject_output = tproject_output)
      DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                       const.dir = f.dir, const.rhs = f.rhs, var = var)
      
    } else {
      
      res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
      
      if (res$status == 0) {
        
        delta <- res$objval
        if (orientation == "oo") {
          delta <- 1 / delta
        }
        res <- res$solution
        
        t <- res[1]
        lambda <- res[2 : (ndr + 1)] / t
        names(lambda) <- dmunames[dmu_ref]
        
        project_input <- res[(ndr + 2) : (ndr + ni + 1)] / t
        names(project_input) <- inputnames
        project_output <- res[(ndr + ni + 2) : (ndr + ni + no + 1)] / t
        names(project_output) <- outputnames
        
        if (compute_target || compute_rho) {
          target_input <- as.vector(inputref %*% lambda)
          names(target_input) <- inputnames
          target_output <- as.vector(outputref %*% lambda)
          names(target_output) <- outputnames
          
          t_input <- project_input - input[, ii] # input superslack
          names(t_input) <- inputnames
          t_output <- output[, ii] - project_output # output superslack
          names(t_output) <- outputnames
          
          slack_input <- project_input - target_input
          names(slack_input) <- inputnames
          slack_output <- target_output - project_output
          names(slack_output) <- outputnames
        }
        
        if (compute_rho) {
          
          input2 <- input
          input2[, ii] <- project_input
          output2 <- output
          output2[, ii] <- project_output
          datadea2 <- structure(list(input = input2,
                                     output = output2,
                                     dmunames = dmunames,
                                     nc_inputs = nc_inputs,
                                     nc_outputs = nc_outputs),
                                class = "deadata")
          
          if ((rts == "grs") && kaizen) {
            stop("Kaizen is not available for generalized returns to scale.")
          }
          
          deasol <- model_sbmeff(datadea = datadea2,
                                 dmu_eval = ii,
                                 dmu_ref = dmu_ref,
                                 weight_input = weight_input[, i],
                                 weight_output = weight_output[, i],
                                 orientation = orientation,
                                 rts = rts,
                                 L = L,
                                 U = U,
                                 kaizen = kaizen,
                                 silent = silent,
                                 compute_target = FALSE)
          
          rho <- deasol$DMU[[1]]$efficiency
          gamma <- rho * delta
          s_input <- deasol$DMU[[1]]$slack_input
          s_output <- deasol$DMU[[1]]$slack_output
          csbm <- (1 - sum((weight_input[, i] * (s_input - t_input)) / (sumwi[i] * input[, ii]))) /
            (1 + sum((weight_output[, i] * (s_output - t_output)) / (sumwo[i] * output[, ii])))
          
        }
        
      } else {
        
        delta <- NA
        lambda <- NA
        project_input <- NA
        project_output <- NA
        if (compute_target || compute_rho) {
          target_input <- NA
          target_output <- NA
          t_input <- NA
          t_output <- NA
          slack_input <- NA
          slack_output <- NA
        }
        if (compute_rho) {
          rho <- NA
          gamma <- NA
          csbm <- NA
        }
        
      }
      
      DMU[[i]] <- list(delta = delta,
                       rho = rho,
                       gamma = gamma,
                       csbm = csbm,
                       kaizen = kaizen,
                       lambda = lambda,
                       project_input = project_input, project_output = project_output,
                       target_input = target_input, target_output = target_output,
                       t_input = t_input, t_output = t_output,
                       slack_input = slack_input, slack_output = slack_output)
      
    }
    
  }
  
  deaOutput <- list(modelname = "sbmsupereff",
                    orientation = orientation,
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    weight_input = weight_input,
                    weight_output = weight_output)
  
  return(structure(deaOutput, class = "dea"))
  
  }
