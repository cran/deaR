#' @title Basic (radial and directional) DEA model.
#'   
#' @description It solves input and output oriented, along with directional, basic
#' DEA models (envelopment form) under constant (CCR model), variable (BCC model),
#' non-increasing, non-decreasing or generalized returns to scale. By default,
#' models are solved in a two-stage process (slacks are maximized).  
#'
#' You can use the \code{model_basic} function to solve directional DEA
#' models by choosing \code{orientation} = "dir".
#' 
#' The model_basic function allows to treat with non-discretional, non-controllable
#' and undesirable inputs/outputs.
#'
#' @note (1) Model proposed by Seiford and Zhu (2002) is applied for undesirable
#' inputs/outputs and non-directional orientation (i.e., input or output oriented).
#' You should select "vrs" returns to scale (BCC model) in order to maintain translation
#' invariance. If deaR detects that you are not specifying \code{rts} = "vrs", it
#' makes the change to "vrs" automatically. 
#' 
#' (2) With undesirable inputs and non-directional orientation use input-oriented
#' BCC model, and with undesirable outputs and non-directional orientation use output-oriented
#' BCC model. Alternatively, you can also treat the undesirable outputs as inputs
#' and then apply the input-oriented BCC model (similarly with undesirable inputs).
#' 
#' (3) Model proposed by Fare and Grosskopf (2004) is applied for undesirable inputs/outputs
#' and directional orientation.
#' 
#' (4) With \code{orientation} = "dir" (directional distance function model), efficient
#' DMUs are those for which \code{beta} = 0.
#' 
#' @usage model_basic(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             orientation = c("io", "oo", "dir"),
#'             dir_input = NULL,
#'             dir_output = NULL,
#'             rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'             L = 1,
#'             U = 1,
#'             maxslack = TRUE,
#'             weight_slack_i = 1,
#'             weight_slack_o = 1,
#'             vtrans_i = NULL,
#'             vtrans_o = NULL,
#'             compute_target = TRUE,
#'             compute_multiplier = FALSE,
#'             returnlp = FALSE,
#'             ...)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s}
#' outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation
#' reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation A string, equal to "io" (input oriented), "oo" (output
#' oriented), or "dir" (directional).
#' @param dir_input A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the input directions.
#' If \code{dir_input} == input matrix (of DMUS in \code{dmu_eval}) and
#' \code{dir_output} == 0, it is equivalent to input oriented (\code{beta} = 1 -
#' \code{efficiency}). If \code{dir_input} is omitted, input matrix (of DMUS in
#' \code{dmu_eval}) is assigned.
#' @param dir_output A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the output directions.
#' If \code{dir_input} == 0 and \code{dir_output} == output matrix (of DMUS in
#' \code{dmu_eval}), it is equivalent to output oriented (\code{beta} = \code{efficiency} - 1).
#' If \code{dir_output} is omitted, output matrix (of DMUS in \code{dmu_eval}) is assigned.
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
#' @param vtrans_i Numeric vector of translation for undesirable inputs with non-directional
#' orientation. If \code{vtrans_i[i]} is \code{NA}, then it applies the "max + 1" translation
#' to the i-th undesirable input. If \code{vtrans_i} is a constant, then it applies
#' the same translation to all undesirable inputs. If \code{vtrans_i} is \code{NULL},
#' then it applies the "max + 1" translation to all undesirable inputs.
#' @param vtrans_o Numeric vector of translation for undesirable outputs with
#' non-directional orientation, analogous to \code{vtrans_i}, but applied to outputs.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets of the
#' max slack solution. 
#' @param compute_multiplier Logical. If it is \code{TRUE}, it computes multipliers
#' (dual solution) when \code{orientation} is "io" or "oo".
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems
#' (objective function and constraints) of stage 1.
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
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1978). “Measuring the efficiency of decision
#' making units”, European Journal of Operational Research 2, 429–444.
#' 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1979). “Short communication: Measuring the
#' efficiency of decision making units”, European Journal of Operational Research 3, 339. 
#' 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1981). "Evaluating Program and Managerial
#' Efficiency: An Application of Data Envelopment Analysis to Program Follow Through",
#' Management Science, 27(6), 668-697. 
#' 
#' Banker, R.; Charnes, A.; Cooper, W.W. (1984). “Some Models for Estimating Technical
#' and Scale Inefficiencies in Data Envelopment Analysis”, Management Science; 30; 1078-1092. 
#'
#' Undesirable inputs/outputs: 
#' 
#' Pastor, J.T. (1996). "Translation Invariance in Data Envelopment Analysis: a
#' Generalization", Annals of Operations Research, 66(2), 91-102. 
#'
#' Seiford, L.M.; Zhu, J. (2002). “Modeling undesirable factors in efficiency evaluation”,
#' European Journal of Operational Research 142, 16-20.
#' 
#' Färe, R. ; Grosskopf, S. (2004). “Modeling undesirable factors in efficiency
#' evaluation: Comment”, European Journal of Operational Research 157, 242-245.
#' 
#' Hua Z.; Bian Y. (2007). DEA with Undesirable Factors. In: Zhu J., Cook W.D. (eds)
#' Modeling Data Irregularities and Structural Complexities in Data Envelopment Analysis.
#' Springer, Boston, MA. 
#' 
#' Non-discretionary/Non-controllable inputs/outputs:  
#'   
#' Banker, R.; Morey, R. (1986). “Efficiency Analysis for Exogenously Fixed Inputs
#' and Outputs”, Operations Research; 34; 513-521. 
#' 
#' Ruggiero J. (2007). Non-Discretionary Inputs. In: Zhu J., Cook W.D. (eds) Modeling
#' Data Irregularities and Structural Complexities in Data Envelopment Analysis.
#' Springer, Boston, MA.
#'
#' Directional DEA model:  
#' 
#' Chambers, R.G.; Chung, Y.; Färe, R. (1996). "Benefit and Distance Functions",
#' Journal of Economic Theory, 70(2), 407-419. 
#' 
#' Chambers, R.G.; Chung, Y.; Färe, R. (1998). "Profit Directional Distance
#' Functions and Nerlovian Efficiency", Journal of Optimization Theory and
#' Applications, 95, 351-354.
#' 
#' @examples 
#' # Example 1. Basic DEA model with desirable inputs/outputs.
#' # Replication of results in Charnes, Cooper and Rhodes (1981).
#' data("PFT1981") 
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ] 
#' PFT <- read_data(PFT, 
#'                  inputs = 2:6, 
#'                  outputs = 7:9 )
#' eval_pft <- model_basic(PFT, 
#'                         orientation = "io", 
#'                         rts = "crs")
#' eff <- efficiencies(eval_pft)
#' s <- slacks(eval_pft) 
#' lamb <- lambdas(eval_pft)
#' tar <- targets(eval_pft)
#' ref <- references(eval_pft) 
#' returns <- rts(eval_pft)
#' 
#' # Example 2. Basic DEA model with undesirable outputs.
#' # Replication of results in Hua and Bian (2007).
#' data("Hua_Bian_2007")
#' # The third output is an undesirable output.
#' data_example <- read_data(Hua_Bian_2007, 
#'                           ni = 2,
#'                           no = 3, 
#'                           ud_outputs = 3) 
#' # Translation parameter (vtrans_o) is set to 1500                          
#' result <- model_basic(data_example, 
#'                       orientation = "oo", 
#'                       rts = "vrs", 
#'                       vtrans_o = 1500) 
#' eff <- efficiencies(result)
#' 1 / eff # results M5 in Table 6-5 (p.119)
#' 
#' # Example 3. Basic DEA model with non-discretionary (fixed) inputs.
#' # Replication of results in Ruggiero (2007).
#' data("Ruggiero2007") 
#' # The second input is a non-discretionary input.
#' datadea <- read_data(Ruggiero2007, 
#'                      ni = 2,
#'                      no = 1, 
#'                      nd_inputs = 2) 
#' result <- model_basic(datadea,
#'                       orientation = "io", 
#'                       rts = "crs")
#' efficiencies(result)
#'  
#' @seealso \code{\link{model_multiplier}}, \code{\link{model_supereff}}
#' 
#' @import lpSolve
#' 
#' @export
  
model_basic <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = c("io", "oo", "dir"),
           dir_input = NULL,
           dir_output = NULL,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           maxslack = TRUE,
           weight_slack_i = 1,
           weight_slack_o = 1,
           vtrans_i = NULL,
           vtrans_o = NULL,
           compute_target = TRUE,
           compute_multiplier = FALSE,
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
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    if (orientation != "dir") {
      datadea_old <- datadea
      res_und <- undesirable_basic(datadea = datadea, vtrans_i = vtrans_i,
                                   vtrans_o = vtrans_o)
      datadea <- res_und$u_datadea
      if (orientation == "oo") {
        vtrans_i <- res_und$vtrans_o
        vtrans_o <- res_und$vtrans_i
      } else {
        vtrans_i <- res_und$vtrans_i
        vtrans_o <- res_und$vtrans_o
      }
      if (!is.null(datadea$ud_inputs) && (orientation != "io")) {
        warning("Undesirable (good) inputs with no input-oriented model.")
      }
      if (!is.null(datadea$ud_outputs) && (orientation != "oo")) {
        warning("Undesirable (bad) outputs with no output-oriented model.")
      }
      if (rts != "vrs") {
        #rts <- "vrs"
        warning("Returns to scale may be changed to variable (vrs) because there
                are data with undesirable inputs/outputs.")
      }
    }
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
  
  if (orientation == "dir") {
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
    if (is.null(dir_input)) {
      dir_input <- matrix(input[, dmu_eval], nrow = ni) # input of DMUs in dmu_eval
    } else {
      if (is.matrix(dir_input)) {
        if ((nrow(dir_input) != ni) || (ncol(dir_input) != nde)) {
          stop("Invalid input direction matrix (number of inputs x number of evaluated DMUs).")
        }
      } else if ((length(dir_input) == 1) || (length(dir_input) == ni)) {
        dir_input <- matrix(dir_input, nrow = ni, ncol = nde)
      } else {
        stop("Invalid input direction vector.")
      }
    }
    rownames(dir_input) <- inputnames
    colnames(dir_input) <- dmunames[dmu_eval]
    
    if (is.null(dir_output)) {
      dir_output <- matrix(output[, dmu_eval], nrow = no) # output of DMUs in dmu_eval
    } else {
      if (is.matrix(dir_output)) {
        if ((nrow(dir_output) != no) || (ncol(dir_output) != nde)) {
          stop("Invalid output direction matrix (number of outputs x number of evaluated DMUs).")
        }
      } else if ((length(dir_output) == 1) || (length(dir_output) == no)) {
        dir_output <- matrix(dir_output, nrow = no, ncol = nde)
      } else {
        stop("Invalid output direction vector.")
      }
    }
    rownames(dir_output) <- outputnames
    colnames(dir_output) <- dmunames[dmu_eval]
    
  } else {
    if (orientation == "io") {
      input <- datadea$input
      output <- datadea$output
      nc_inputs <- datadea$nc_inputs
      nc_outputs <- datadea$nc_outputs
      nd_inputs <- datadea$nd_inputs
      nd_outputs <- datadea$nd_outputs
      ud_inputs <- datadea$ud_inputs
      ud_outputs <- datadea$ud_outputs
      obj <- "min"
      orient <- 1
    } else {
      input <- -datadea$output
      output <- -datadea$input
      nc_inputs <- datadea$nc_outputs
      nc_outputs <- datadea$nc_inputs
      nd_inputs <- datadea$nd_outputs
      nd_outputs <- datadea$nd_inputs
      ud_inputs <- datadea$ud_outputs
      ud_outputs <- datadea$ud_inputs
      aux <- weight_slack_i
      weight_slack_i <- weight_slack_o
      weight_slack_o <- aux
      obj <- "max"
      orient <- -1
    }
    inputnames <- rownames(input)
    outputnames <- rownames(output)
    ni <- nrow(input) # number of  inputs
    no <- nrow(output) # number of outputs
  }
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  ncd_inputs <- c(nc_inputs, nd_inputs)
  
  target_input <- NULL
  target_output <- NULL
  multiplier_input <- NULL
  multiplier_output <- NULL
  multiplier_rts <- NULL
  orientation_param <- NULL
  
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
    
    ### If orientation == "dir", slacks for undesirable i/o are 0 ###
    if (orientation == "dir") {
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
  }
  
  if (orientation != "dir") { ############## orientation != "dir" ##############
    
    # Constraints matrix of 2nd bloc of constraints stage 1
    f.con.2 <- cbind(matrix(0, nrow = no, ncol = 1), outputref)
    
    for (i in 1:nde) {
      
      ii <- dmu_eval[i]
      
      # Constraints matrix stage 1
      f.con.1 <- cbind(-input[, ii], inputref)
      f.con.1[ncd_inputs, 1] <- 0
      f.con <- rbind(f.con.1, f.con.2, f.con.rs)
      
      # Right hand side vector stage 1
      f.rhs <- c(rep(0, ni), output[, ii], f.rhs.rs)
      f.rhs[ncd_inputs] <- input[ncd_inputs, ii]
      
      if (returnlp) {
        
        lambda <- rep(0, ndr)
        names(lambda) <- dmunames[dmu_ref]
        var <- list(efficiency = 0, lambda = lambda)
        DMU[[i]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                         const.dir = f.dir, const.rhs = f.rhs, var = var)
        
      } else {
        
        if (compute_multiplier) {
          
          res <- lp(obj, f.obj, f.con, f.dir, f.rhs, compute.sens = TRUE)
          
          if (res$status == 0) {
            multiplier_input <- -orient * res$duals[1 : ni]
            names(multiplier_input) <- inputnames
            multiplier_output <- orient * res$duals[(ni + 1) : (ni + no)]
            names(multiplier_output) <- outputnames
            if (rts == "grs") {
              multiplier_rts <- res$duals[(ni + no + 1):(ni + no + 2)]
              names(multiplier_rts) <- c("grs_L", "grs_U")
            } else if (rts != "crs") {
              multiplier_rts <- res$duals[ni + no + 1]
              names(multiplier_rts) <- rts
            }
          } else {
            multiplier_input <- NA
            multiplier_output <- NA
            multiplier_rts <- NA
          }
          
        } else {
          res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
        }
        
        if(res$status == 0) {
          
          res <- res$solution
          
          efficiency <- res[1]
          
          if (maxslack) {
            
            # Objective function coefficients stage 2
            f.obj2 <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i])
            
            # Right hand side vector stage 2
            f.rhs2 <- c(efficiency * input[, ii], output[, ii], rep(0, nnci + nnco), f.rhs.rs)
            f.rhs2[ncd_inputs] <- input[ncd_inputs, ii]
            
            res <- lp("max", f.obj2, f.con2, f.dir2, f.rhs2)$solution
            
            lambda <- res[1 : ndr]
            names(lambda) <- dmunames[dmu_ref]
            
            slack_input <- res[(ndr + 1) : (ndr + ni)]
            names(slack_input) <- inputnames
            slack_output <- res[(ndr + ni + 1) : (ndr + ni + no)]
            names(slack_output) <- outputnames
            
            if (compute_target) {
              target_input <- orient * as.vector(inputref %*% lambda)
              target_output <- orient * as.vector(outputref %*% lambda)
              #target_input <- orient * (efficiency * input[, ii] - slack_input) # Alternative
              names(target_input) <- inputnames
              #target_output <- orient * (output[, ii] + slack_output) # Alternative
              names(target_output) <- outputnames
              target_input[ud_inputs] <- vtrans_i - target_input[ud_inputs]
              target_output[ud_outputs] <- vtrans_o - target_output[ud_outputs]
            }
            
          } else {
            
            lambda <- res[2 : (ndr + 1)]
            names(lambda) <- dmunames[dmu_ref]
            
            target_input <- orient * as.vector(inputref %*% lambda)
            names(target_input) <- inputnames
            target_output <- orient * as.vector(outputref %*% lambda)
            names(target_output) <- outputnames
            
            slack_input <- efficiency * input[, ii] - orient * target_input
            slack_input[ncd_inputs] <- input[ncd_inputs, ii] - orient * target_input[ncd_inputs]
            names(slack_input) <- inputnames
            slack_output <- orient * target_output - output[, ii]
            names(slack_output) <- outputnames
            
            target_input[ud_inputs] <- vtrans_i - target_input[ud_inputs]
            target_output[ud_outputs] <- vtrans_o - target_output[ud_outputs]
            
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
        
        if (orientation == "io") {
          if (rts == "crs") {
            DMU[[i]] <- list(efficiency = efficiency,
                             lambda = lambda,
                             slack_input = slack_input, slack_output = slack_output,
                             target_input = target_input, target_output = target_output,
                             multiplier_input = multiplier_input, multiplier_output = multiplier_output)
          } else {
            DMU[[i]] <- list(efficiency = efficiency,
                             lambda = lambda,
                             slack_input = slack_input, slack_output = slack_output,
                             target_input = target_input, target_output = target_output,
                             multiplier_input = multiplier_input, multiplier_output = multiplier_output,
                             multiplier_rts = multiplier_rts)
          }
        } else {
          aux <- weight_slack_i
          weight_slack_i <- weight_slack_o
          weight_slack_o <- aux
          if (rts == "crs") {
            DMU[[i]] <- list(efficiency = efficiency, # 1/ efficiency alternative
                             lambda = lambda,
                             slack_input = slack_output, slack_output = slack_input,
                             target_input = target_output, target_output = target_input,
                             multiplier_input = multiplier_output, multiplier_output = multiplier_input)
          } else {
            DMU[[i]] <- list(efficiency = efficiency, # 1/ efficiency alternative
                             lambda = lambda,
                             slack_input = slack_output, slack_output = slack_input,
                             target_input = target_output, target_output = target_input,
                             multiplier_input = multiplier_output, multiplier_output = multiplier_input,
                             multiplier_rts = multiplier_rts)
          }
        }
        
      }
      
    }
    
  } else {                    ############## orientation == "dir" ##############
    
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
            
            if (compute_target) {
              target_input <- as.vector(inputref %*% lambda)
              target_output <- as.vector(outputref %*% lambda)
              names(target_input) <- inputnames
              names(target_output) <- outputnames
            }
            
          } else {
            
            lambda <- res[2 : (ndr + 1)]
            names(lambda) <- dmunames[dmu_ref]
            
            target_input <- as.vector(inputref %*% lambda)
            names(target_input) <- inputnames
            target_output <- as.vector(outputref %*% lambda)
            names(target_output) <- outputnames
            
            slack_input <- input[, ii] - beta * dir_input[, i] - target_input
            slack_input[ud_inputs] <- input[ud_inputs, ii] + beta * dir_input[ud_inputs, i] -
              target_input[ud_inputs]
            slack_input[ncd_inputs] <- input[ncd_inputs, ii] - target_input[ncd_inputs]
            names(slack_input) <- inputnames
            slack_output <- target_output - output[, ii] - beta * dir_output[, i]
            slack_output[ud_outputs] <- target_output[ud_outputs] - output[ud_outputs, ii] +
              beta * dir_output[ud_outputs, i]
            slack_output[ncd_outputs] <- target_output[ncd_outputs] - output[ncd_outputs, ii]
            names(slack_output) <- outputnames
            
          }
          
        } else {
          
          beta <- NA
          lambda <- NA
          slack_input <- NA
          slack_output <- NA
          if (compute_target) {
            target_input <- NA
            target_output <- NA
          }
          
        }
        
        DMU[[i]] <- list(beta = beta,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output
                         )
        
      }
      
    }
    
    orientation_param <- list(
                        dir_input = dir_input,
                        dir_output = dir_output)
    
  }
  
  if ((!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) && (orientation != "dir")) {
    datadea <- datadea_old
    vtrans_i <- res_und$vtrans_i
    vtrans_o <- res_und$vtrans_o
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
  
  deaOutput <- list(modelname = "basic",
                   orientation = orientation,
                   orientation_param = orientation_param,
                   rts = rts,
                   L = L,
                   U = U,
                   DMU = DMU,
                   data = datadea,
                   dmu_eval = dmu_eval,
                   dmu_ref = dmu_ref,
                   vtrans_i = vtrans_i,
                   vtrans_o = vtrans_o,
                   maxslack = maxslack,
                   weight_slack_i = weight_slack_i,
                   weight_slack_o = weight_slack_o)
 
  return(structure(deaOutput, class = "dea"))
}
