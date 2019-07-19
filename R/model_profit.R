#' @title Profit efficiency DEA model.
#'   
#' @description Cost, revenue and profit efficiency DEA models.
#' @usage model_profit(datadea,
#'              dmu_eval = NULL,
#'              dmu_ref = NULL,
#'              price_input = NULL,
#'              price_output = NULL,
#'              rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'              L = 1,
#'              U = 1,
#'              restricted_optimal = TRUE,
#'              returnlp = FALSE,
#'              ...)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' @param price_input Unit prices of inputs for cost or profit efficiency models. It is a value, vector of length \code{m}, or matrix \code{m} x \code{ne} (where \code{ne} is the lenght of \code{dmu_eval}).
#' @param price_output Unit prices of outputs for revenue or profit efficiency models. It is a value, vector of length \code{s}, or matrix \code{s} x \code{ne}.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#'            "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param restricted_optimal Logical. If it is \code{TRUE}, the optimal inputs are restricted to be <= inputs (for cost efficiency models) or
#'                           the optimal outputs are restricted to be >= outputs (for revenue efficiency models).
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective function and constraints) of stage 1.
#' @param ... Ignored, for compatibility issues.
#'
#' @references
#' Coelli, T.; Prasada Rao, D.S.; Battese, G.E. (1998). An introduction to efficiency and productivity analysis. Jossey-Bass, San Francisco, pp 73–104. \url{https://doi.org/10.1002/ev.1441}  
#'    
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
#' Coelli, T.; Prasada Rao, D.S.; Battese, G.E. An introduction to efficiency and productivity analysis. Boston: Kluwer Academic Publishers.
#' 
#' @examples 
#' # Example 1. Replication of results in Coelli et al. (1998, p.166).
#' # Cost efficiency model.
#' data("Coelli_1998")
#' # Selection of prices: input_prices is the transpose where the prices for inputs are. 
#' input_prices <- t(Coelli_1998[, 5:6]) 
#' 
#' data_example1 <- read_data(Coelli_1998,
#'                            ni = 2,
#'                            no = 1)
#' result1 <- model_profit(data_example1,
#'                        price_input = input_prices,
#'                        rts = "crs", 
#'                        restricted_optimal = FALSE) 
#' # notice that the option by default is restricted_optimal = TRUE
#' efficiencies(result1)
#' 
#' # Example 2. Revenue efficiency model.
#' data("Coelli_1998")
#' # Selection of prices for output: output_prices is the transpose where the prices for outputs are. 
#' output_prices <- t(Coelli_1998[, 7]) 
#' data_example2 <- read_data(Coelli_1998,
#'                           ni = 2,
#'                           no = 1)
#' result2 <- model_profit(data_example2,
#'                        price_output = output_prices,
#'                        rts = "crs", 
#'                        restricted_optimal = FALSE) 
#' # notice that the option by default is restricted_optimal = TRUE
#' efficiencies(result2)
#' 
#' # Example 3. Profit efficiency model.
#' data("Coelli_1998")
#' # Selection of prices for inputs and outputs: input_prices and output_prices are 
#' # the transpose where the prices (for inputs and outputs) are. 
#' input_prices <- t(Coelli_1998[, 5:6]) 
#' output_prices <- t(Coelli_1998[, 7]) 
#' data_example3 <- read_data(Coelli_1998,
#'                            ni = 2,
#'                            no = 1)
#' result3 <- model_profit(data_example3,
#'                         price_input = input_prices,
#'                         price_output = output_prices,
#'                         rts = "crs", 
#'                         restricted_optimal = FALSE) 
#' # notice that the option by default is restricted_optimal = TRUE
#' efficiencies(result3)
#' 
#' @seealso \code{\link{model_deaps}}, \code{\link{model_nonradial}}, \code{\link{model_sbmeff}}
#'  
#' @import lpSolve
#' 
#' @export
  
model_profit <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           price_input = NULL,
           price_output = NULL,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,           
           U = 1,
           restricted_optimal = TRUE,
           returnlp = FALSE,
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
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
  nc_inputs <- datadea$nc_inputs
  nc_outputs <- datadea$nc_outputs
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  aux_i <- 1
  aux_o <- 1
  # Checking prices
  if (is.null(price_input)) {
    aux_i = 0
  } else {
    if (is.matrix(price_input)) {
      if ((nrow(price_input) != ni) || (ncol(price_input) != nde)) {
        stop("Invalid input prices matrix (number of inputs x number of evaluated DMUs).")
      }
    } else if ((length(price_input) == 1) || (length(price_input) == ni)) {
      price_input <- matrix(price_input, nrow = ni, ncol = nde)
    } else {
      stop("Invalid input prices vector (number of inputs).")
    }
    price_input[nc_inputs, ] <- 0
    rownames(price_input) <- inputnames
    colnames(price_input) <- dmunames[dmu_eval]
  }
  if (is.null(price_output)) {
    aux_o = 0
  } else {
    if (is.matrix(price_output)) {
      if ((nrow(price_output) != no) || (ncol(price_output) != nde)) {
        stop("Invalid output prices matrix (number of outputs x number of evaluated DMUs).")
      }
    } else if ((length(price_output) == 1) || (length(price_output) == no)) {
      price_output <- matrix(price_output, nrow = no, ncol = nde)
    } else {
      stop("Invalid output prices vector (number of outputs).")
    }
    price_output[nc_outputs, ] <- 0
    rownames(price_output) <- outputnames
    colnames(price_output) <- dmunames[dmu_eval]
  }
  if ((aux_i == 0) && (aux_o == 0)) {
    stop("No prices specified.")
  }
  
  ### restricted_optimal must be TRUE for profit efficiency models ###
  if (!restricted_optimal && (aux_i * aux_o == 1)) {
    restricted_optimal <- TRUE
  }
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  ###########################
  
  if (rts == "crs") {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- cbind(matrix(0, nrow = 1, ncol = ni + no), matrix(1, nrow = 1, ncol = ndr))
    f.rhs.rs <- 1
    if (rts == "vrs") {
      f.dir.rs <- "="
    } else if (rts == "nirs") {
      f.dir.rs <- "<="
    } else if (rts == "ndrs") {
      f.dir.rs <- ">="
    } else {
      f.con.rs <- rbind(f.con.rs, f.con.rs)
      f.dir.rs <- c(">=", "<=")
      f.rhs.rs <- c(L, U)
    }
  }
  
  # Matriz técnica
  f.con.1 <- cbind(-diag(ni), matrix(0, nrow = ni, ncol = no), inputref)
  f.con.2 <- cbind(matrix(0, nrow = no, ncol = ni), -diag(no), outputref)
  f.con.3 <- cbind(diag(ni), matrix(0, nrow = ni, ncol = no + ndr))
  f.con.4 <- cbind(matrix(0, nrow = no, ncol = ni), diag(no), matrix(0, nrow = no, ncol = ndr))
  if (!restricted_optimal) {
    if (aux_i == 0) {
      f.con.4 <- NULL
    } else if (aux_o == 0) {
      f.con.3 <- NULL
    }
  }
  f.con <- rbind(f.con.1, f.con.2, f.con.3, f.con.4, f.con.rs)
  
  # Vector de dirección de restricciones
  if (aux_i == 0) {
    f.dir.3 <- rep("=", ni)
  } else {
    f.dir.3 <- rep("<=", ni)
    f.dir.3[nc_inputs] <- "="
  }
  if (aux_o == 0) {
    f.dir.4 <- rep("=", no)
  } else {
    f.dir.4 <- rep(">=", no)
    f.dir.4[nc_outputs] <- "="
  }
  if (!restricted_optimal) {
    if (aux_i == 0) {
      f.dir.4 <- NULL
    } else if (aux_o == 0) {
      f.dir.3 <- NULL
    }
  }
  f.dir <- c(rep("<=", ni), rep(">=", no), f.dir.3, f.dir.4, f.dir.rs)
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    # Vector de coeficientes de la función objetivo
    if (aux_o == 0) {
      f.obj <- c(-price_input[, i], rep(0, no + ndr))
    } else if (aux_i == 0) {
      f.obj <- c(rep(0, ni), price_output[, i], rep(0, ndr))
    } else {
      f.obj <- c(-price_input[, i], price_output[, i], rep(0, ndr))
    }
    
    # Vector de términos independientes
    f.rhs.3 <- input[, ii]
    f.rhs.4 <- output[, ii]
    if (!restricted_optimal) {
      if (aux_i == 0) {
        f.rhs.4 <- NULL
      } else if (aux_o == 0) {
        f.rhs.3 <- NULL
      }
    }
    f.rhs <- c(rep(0, ni + no), f.rhs.3, f.rhs.4, f.rhs.rs)
    
    if (returnlp) {
      
      optimal_input = rep(0, ni)
      names(optimal_input) <- inputnames
      optimal_output = rep(0, no)
      names(optimal_output) <- outputnames
      lambda <- rep(0, ndr)
      names(lambda) <- dmunames[dmu_ref]
      var <- list(optimal_input = optimal_input, optimal_output = optimal_output, lambda = lambda)
      DMU[[i]] <- list(direction = "max", objective.in = f.obj, const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs,
                       var = var)
      
    } else {
      
      res <- lp("max", f.obj, f.con, f.dir, f.rhs)
      
      if (res$status == 0) {
        
        objval <- res$objval
        
        res <- res$solution
        
        optimal_input <- res[1 : ni]
        names(optimal_input) <- inputnames
        optimal_output <- res[(ni + 1) : (ni + no)]
        names(optimal_output) <- outputnames
        lambda <- res[(ni + no + 1) : (ni + no + ndr)]
        names(lambda) <- dmunames[dmu_ref]
        
        if (aux_o == 0) {
          eff <- (price_input[, i] %*% optimal_input) / 
            (price_input[, i] %*% input[, ii])
          objval <- -objval
        } else if (aux_i == 0)
          eff <- (price_output[, i] %*% output[, ii]) / 
          (price_output[, i] %*% optimal_output)
        else {
          eff <- (-price_input[, i] %*% input[, ii] + price_output[, i] %*% output[, ii]) / 
            (-price_input[, i] %*% optimal_input + price_output[, i] %*% optimal_output)
        }
        
        
      } else {
        
        objval <- NA
        optimal_input <- NA
        optimal_output <- NA
        lambda <- NA
        eff <- NA
        
      }
      
      DMU[[i]] <- list(objval = objval,
                       efficiency = eff,
                       lambda = lambda,
                       optimal_input = optimal_input, optimal_output = optimal_output)
      
    }
    
  }
 
  deaOutput <- list(modelname = "profit",
                    rts = rts,
                    L = L,
                    U = U,
                    DMU = DMU,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    price_input = price_input,
                    price_output = price_output,
                    restricted_optimal = restricted_optimal,
                    orientation = "io")
 
  return(structure(deaOutput, class = "dea"))
 
}
