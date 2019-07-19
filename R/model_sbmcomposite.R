#' @title SBM composite superefficiency model
#'   
#' @description SBM composite superefficiency model.
#' 
#' @usage model_sbmcomposite(datadea,
#'                    dmu_eval = NULL,
#'                    dmu_ref = NULL,
#'                    weight_input = 1,
#'                    weight_output = 1,
#'                    orientation = c("no", "io", "oo"),
#'                    rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                    L = 1,
#'                    U = 1,
#'                    compute_target = TRUE,
#'                    compute_alt = FALSE,
#'                    delta_margin = 1e-8,
#'                    returnlp = FALSE)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' @param weight_input A value, vector of length \code{m}, or matrix \code{m} x \code{ne} (where \code{ne} is the lenght of \code{dmu_eval})
#' with weights to inputs corresponding to the relative importance of items.
#' @param weight_output A value, vector of length \code{m}, or matrix \code{m} x \code{ne} (where \code{ne} is the lenght of \code{dmu_eval})
#' with weights to outputs corresponding to the relative importance of items.
#' @param orientation A string, equal to "no" (non-oriented), "io" (input-oriented) or "oo" (output-oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets, superslacks (\code{t_input} and \code{t_output}) and slacks.
#' @param compute_alt Logical. If it is \code{TRUE}, it computes the alternative SBM composite super-efficiency.
#' @param delta_margin Some numerical errors may appear due to truncation of the SBM super-efficiency in some DMUs.
#' In these cases, the SBM super-efficiency is slightly increased by multiplying it by 1 + \code{delta_margin}.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective function and constraints).
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
#' Tone, K. (2002). "A slacks-based measure of super-efficiency in data envelopment analysis", European Journal of Operational Research, 143, 32-41. \url{https://doi.org/10.1016/S0377-2217(01)00324-1}  
#' 
#' Cooper, W.W.; Seiford, L.M.; Tone, K. (2007). Data Envelopment Analysis. A Comprehensive Text with Models, Applications, References and DEA-Solver Software. 2nd Edition. Springer, New York. \url{https://doi.org/10.1007/978-0-387-45283-8}
#' 
#' @examples 
#' # Data from Tone(2002, p.39)
#' data("Power_plants")
#' data_example <- read_data(Power_plants,
#'                           ni = 4,
#'                           no = 2)
#' result <- model_sbmcomposite(data_example,
#'                              compute_alt = TRUE) 
#' 
#' @seealso \code{\link{model_sbmsupereff}}, \code{\link{model_sbmeff}}, \code{\link{model_supereff}}, \code{\link{model_addsupereff}}
#' 
#' @import lpSolve
#' 
#' @export

model_sbmcomposite <-
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
           compute_alt = FALSE,
           delta_margin = 1e-8,
           returnlp = FALSE) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
  # Checking non-controllable or non-discretionary inputs/outputs
  #if ((!is.null(datadea$nc_inputs)) || (!is.null(datadea$nc_outputs))
  #    || (!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
  #  warning("This model does not take into account non-controllable or non-discretionary feature for inputs/outputs.")
  #  datadea$nc_inputs <- NULL
  #  datadea$nc_outputs <- NULL
  #  datadea$nd_inputs <- NULL
  #  datadea$nd_outputs <- NULL
  #}
  
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
  
  ########## Compute SBM superefficiency ##########
  
  deasol <- model_sbmsupereff(datadea = datadea,
                              dmu_eval = dmu_eval,
                              dmu_ref = dmu_ref,
                              weight_input = weight_input,
                              weight_output = weight_output,
                              orientation = orientation,
                              rts = rts,
                              L = L,
                              U = U,
                              compute_target = FALSE)
  
  delta <- unlist(lapply(deasol$DMU, function(x)
    x$delta))
  
  if (any(is.na(delta))) {
    stop("There are some errors in the SBM super-efficiency model.")
  }
  
  #################################################
  
  target_input <- NULL
  target_output <- NULL
  project_input <- NULL
  project_output <- NULL
  slack_input <- NULL
  slack_output <- NULL
  superslack_input <- NULL
  superslack_output <- NULL
  gamma <- NULL
  gamma_alt <- NULL
  rho <- NULL
    
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
    
  ###########################
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else if (rts == "grs") {
    f.con.rs <- rbind(cbind(-L, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2 * ni + 2 * no)),
                      cbind(-U, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2 * ni + 2 * no)))
    f.dir.rs <- c(">=", "<=")
    f.rhs.rs <- c(0, 0)
  } else {
    f.con.rs <- cbind(-1, matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = 2 * ni + 2 * no))
    f.rhs.rs <- 0
    if (rts == "vrs") {
      f.dir.rs <- "="
    } else if (rts == "nirs") {
      f.dir.rs <- "<="
    } else {
      f.dir.rs <- ">="
    }
  }
  
  # Vector de dirección de restricciones, matriz técnica y términos independientes
  if (orientation == "no") {
    f.con.or <- NULL
    f.dir <- c(rep("=", 3 + ni + no), f.dir.rs)
    f.rhs <- c(1, rep(0, ni + no + 2), f.rhs.rs)
  } else if (orientation == "io") {
    f.con.0 <- c(1, rep(0, ndr + 2 * ni + 2 * no))
    f.con.or <- cbind(matrix(0, nrow = 2 * no, ncol = 1 + ndr + 2 * ni), diag(2 * no))
    f.dir <- c(rep("=", 2 + ni), rep (">=", no), rep("=", 2 * no + 1), f.dir.rs)
  } else {
    f.con.0 <- c(1, rep(0, ndr + 2 * ni + 2 * no))
    f.con.or <- cbind(matrix(0, nrow = 2 * ni, ncol = 1 + ndr), diag(2 * ni), matrix(0, nrow = 2 * ni, ncol = 2 * no))
    f.dir <- c(rep("=", 2), rep ("<=", ni), rep ("=", 2 * ni + no + 1), f.dir.rs)
  }
  f.dir[2 + c(nc_inputs, ni + nc_outputs)] <- "="
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    # Vector de coeficientes de la función objetivo, Matriz técnica y Vector de dirección de restricciones
    if (orientation == "no") {
      auxi <- weight_input[, i] / (sumwi[i] * input[, ii])
      auxo <- weight_output[, i] / (sumwo[i] * output[, ii])
      f.obj <- c(1, rep(0, ndr), -auxi, auxi, rep(0, 2 * no))
      f.con.0 <- c(1, rep(0, ndr + 2 * ni), auxo, -auxo)
      f.con.00 <- c(1 - delta[i], rep(0, ndr + ni), auxi, rep(0, no), delta[i] * auxo)
      f.con.1 <- cbind(-input[, ii], inputref, diag(ni), -diag(ni), matrix(0, nrow = ni, ncol = 2 * no))
      f.con.2 <- cbind(-output[, ii], outputref, matrix(0, nrow = no, ncol = 2 * ni), -diag(no), diag(no))
    } else if (orientation == "io") {
      auxi <- weight_input[, i] / (sumwi[i] * input[, ii])
      f.obj <- c(1, rep(0, ndr), -auxi, auxi, rep(0, 2 * no))
      f.con.00 <- c(0, rep(0, ndr + ni), auxi, rep(0, 2 * no))
      f.con.1 <- cbind(0, inputref, diag(ni), -diag(ni), matrix(0, nrow = ni, ncol = 2 * no))
      f.con.2 <- cbind(0, outputref, matrix(0, nrow = no, ncol = 2 * (ni + no)))
      f.rhs <- c(1, delta[i] - 1, input[, ii], output[, ii], rep(0, 2 * no + 1), f.rhs.rs)
    } else {
      auxo <- weight_output[, i] / (sumwo[i] * output[, ii])
      f.obj <- c(1, rep(0, ndr + 2 * ni), auxo, -auxo)
      f.con.00 <- c(0, rep(0, ndr + 2 * ni + no), -auxo)
      f.con.1 <- cbind(0, inputref, matrix(0, nrow = ni, ncol = 2 * (ni + no)))
      f.con.2 <- cbind(0, outputref, matrix(0, nrow = no, ncol = 2 * ni), -diag(no), diag(no))
      f.rhs <- c(1, 1 / delta[i] - 1, input[, ii], output[, ii], rep(0, 2 * ni + 1), f.rhs.rs)
    }
    
    # Matriz técnica
    f.con.se <- rep(0, ndr)
    f.con.se[dmu_ref == ii] <- 1
    f.con.se <- c(0, f.con.se, rep(0, 2 * (ni + no)))
    f.con <- rbind(f.con.0, f.con.00, f.con.1, f.con.2, f.con.or, f.con.se, f.con.rs)
    
    if (returnlp) {
      
      t <- 0
      names(t) <- "t"
      tlambda <- rep(0, ndr)
      names(tlambda) <- dmunames[dmu_ref]
      tslack_input <- rep(0, ni)
      names(tslack_input) <- inputnames
      tsuperslack_input <- tslack_input
      tslack_output <- rep(0, no)
      names(tslack_output) <- outputnames
      tsuperslack_output <- tslack_output
      var <- list(t = t, tlambda = tlambda, tslack_input = tslack_input, tsuperslack_input = tsuperslack_input,
                  tslack_output = tslack_output, tsuperslack_output = tsuperslack_output)
      DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs,
                       var = var)
      
    } else {
      
      res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
      
      if (res$status != 0) { # Second opportunity: apply delta_margin to prevent NA due to truncation
        
        delta[i] <- delta[i] * (1 + delta_margin)
        if (orientation == "no") {
          f.con.00 <- c(1 - delta[i], rep(0, ndr + ni), auxi, rep(0, no), delta[i] * auxo)
          f.con <- rbind(f.con.0, f.con.00, f.con.1, f.con.2, f.con.or, f.con.se, f.con.rs)
        } else if (orientation == "io") {
          f.rhs <- c(1, delta[i] - 1, input[, ii], output[, ii], rep(0, 2 * no + 1), f.rhs.rs)
        } else {
          f.rhs <- c(1, 1 / delta[i] - 1, input[, ii], output[, ii], rep(0, 2 * ni + 1), f.rhs.rs)
        }
        res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
        
      }
      
      if (res$status == 0) {
        
        gamma <- res$objval
        if (orientation == "oo") {
          gamma <- 1 / gamma
        }
        res <- res$solution
        
        t <- res[1]
        lambda <- res[2 : (ndr + 1)] / t
        names(lambda) <- dmunames[dmu_ref]
        
        slack_input <- res[(ndr + 2) : (ndr + ni + 1)] / t
        names(slack_input) <- inputnames
        superslack_input <- res[(ndr + ni + 2) : (ndr + 2 * ni + 1)] / t
        names(superslack_input) <- inputnames
        slack_output <- res[(ndr + 2 * ni + 2) : (ndr + 2 * ni + no + 1)] / t
        names(slack_output) <- outputnames
        superslack_output <- res[(ndr + 2 * ni + no + 2) : (ndr + 2 * ni + 2 * no + 1)] / t
        names(superslack_output) <- outputnames
        
        project_input <- input[, ii] + superslack_input
        project_output <- output[, ii] - superslack_output
        
        if (compute_target) {
          target_input <- as.vector(inputref %*% lambda)
          names(target_input) <- inputnames
          target_output <- as.vector(outputref %*% lambda)
          names(target_output) <- outputnames
        }
        
        if (compute_alt) {
          
          if (orientation == "no") {
            f.obj <- c(1, rep(0, ndr), -auxi, rep(0, ni + 2 * no))
            f.con.0 <- c(1, rep(0, ndr + 2 * ni), auxo, rep(0, no))
            f.con <- rbind(f.con.0, f.con.00, f.con.1, f.con.2, f.con.or, f.con.se, f.con.rs) 
          } else if (orientation == "io") {
            f.obj <- c(1, rep(0, ndr), -auxi, rep(0, ni + 2 * no))
          } else {
            f.obj <- c(1, rep(0, ndr + 2 * ni), auxo, rep(0, no))
          }
          
          res_alt <- lp(obj, f.obj, f.con, f.dir, f.rhs)
          rho <- res_alt$objval
          if (orientation == "oo") {
            rho <- 1 / rho
          }
          gamma_alt <- delta[i] + rho - 1
          
        }
        
      } else {
        
        warning(c("There are numerical errors due to truncation in DMU ",
                dmunames[i], ", try to increase delta_margin."))
        
        gamma <- NA
        lambda <- NA
        slack_input <- NA
        superslack_input <- NA
        slack_output <- NA
        superslack_output <- NA
        project_input <- NA
        project_output <- NA
        if (compute_target) {
          target_input <- NA
          target_output <- NA
        }
        if (compute_alt) {
          gamma_alt <- NA
          rho <- NA
        }
        
      }
      
      DMU[[i]] <- list(gamma = gamma,
                       delta = delta[i],
                       rho = rho,
                       gamma_alt = gamma_alt,
                       t = t,
                       lambda = lambda,
                       slack_input = slack_input, superslack_input = superslack_input,
                       slack_output = slack_output, superslack_output = superslack_output,
                       project_input = project_input, project_output = project_output,
                       target_input = target_input, target_output = target_output)
      
    }
    
  }
  
  deaOutput <- list(modelname = "sbmcomposite",
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
