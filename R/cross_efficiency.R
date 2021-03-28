#' @title Cross efficiency tables
#'   
#' @description Computes arbitrary, benevolent and aggressive formulation of cross-efficiency under constant and variable returns-to-scale. Doyle and Green (1994) present three alternatives ways of formulating the secondary goal (wich will minimize or maximize the other DMUs' cross-efficiencies in some way). Methods II and III are implemented in deaR with constant returns-to-scale. The maverick index is also calculated. 
#' @note (1) We can obtain negative cross-efficiency in the input-oriented DEA model under variable returns to scale. However, the same does not happen in the case of the output oriented VRS DEA model. For this reason, the proposal of Lim and Zhu (2015) is implemented in deaR to calculate the input-oriented cross-efficiency model under variable returns-to-scale.
#' 
#' (2) The multiplier model can have alternate optimal solutions (see note 1 in model_multiplier). So, depending on the optimal weights selected we can obtain different cross-efficinecy scores.
#' 
#' @usage cross_efficiency(datadea,
#'                  dmu_eval = NULL,
#'                  dmu_ref = NULL,
#'                  epsilon = 0, 
#'                  orientation = c("io", "oo"),
#'                  rts = c("crs", "vrs"),
#'                  selfapp = TRUE,
#'                  correction = FALSE,
#'                  M2 = TRUE,
#'                  M3 = TRUE)
#' 
#' @param datadea An object of class \code{dea} or \code{deadata}. If it is of class \code{dea} it must have
#'                been obtained with some of the multiplier DEA models.
#' @param dmu_eval A numeric vector. Only the multipliers of DMUs in \code{dmu_eval} are computed.
#'                 If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#'                If \code{NULL} (default), all DMUs are considered.
#' @param epsilon Numeric, multipliers must be >= \code{epsilon}.
#' @param orientation A string, equal to "io" (input-oriented) or "oo" (output-oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant) or
#'            "vrs" (variable).
#' @param selfapp Logical. If it is \code{TRUE}, self-appraisal is included in the average scores of
#'                \code{A} and \code{e}.
#' @param correction Logical. If it is \code{TRUE}, a correction is applied in the "vrs" input-oriented
#' model in order to avoid negative cross-efficiencies, according to Lim & Zhu (2015).
#' @param M2 Logical. If it is \code{TRUE}, it computes Method II for aggresive/benevolent estimations.
#' @param M3 Logical. If it is \code{TRUE}, it computes Method III for aggresive/benevolent estimations.
#'   
#'
#' @references
#' Sexton, T.R., Silkman, R.H.; Hogan, A.J. (1986). Data envelopment analysis: critique and extensions. In: Silkman RH (ed) Measuring efficiency: an assessment of data envelopment analysis, vol 32. Jossey-Bass, San Francisco, pp 73–104. \doi{10.1002/ev.1441}  
#' 
#' Doyle, J.; Green, R. (1994). “Efficiency and cross efficiency in DEA: derivations, meanings and the uses”,  Journal of Operational Research Society, 45(5), 567–578. DOI: 10.2307/2584392 
#'  
#' Cook, W.D.; Zhu, J. (2015). DEA Cross Efficiency. In: Zhu, J. (ed) Data Envelopment Analysis. A Handbook of Models and Methods. International Series in Operations Research & Management Science, vol 221. Springer, Boston, MA, 23-43 \doi{10.1007/978-1-4899-7553-9_2} 
#'  
#' Lim, S.; Zhu, J. (2015). "DEA Cross-Efficiency Under Variable Returns to Scale". Journal of Operational Research Society, 66(3), p. 476-487. \doi{10.1057/jors.2014.13}
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
#' @examples 
#' # Example 1.
#' # Arbitrary formulation. Input-oriented model under constant returns-to-scale.
#' data("Golany_Roll_1989")
#' data_example <- read_data(datadea = Golany_Roll_1989, 
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- cross_efficiency(data_example, 
#'                            orientation = "io", 
#'                            rts = "crs", 
#'                            selfapp = TRUE)
#' result$Arbitrary$cross_eff
#' result$Arbitrary$e
#' 
#' # Example 2.
#' # Benevolent formulation (method II). Input-oriented.
#' data("Golany_Roll_1989")
#' data_example <- read_data(datadea = Golany_Roll_1989, 
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- cross_efficiency(data_example, 
#'                            orientation = "io", 
#'                            selfapp = TRUE)
#' result$M2_ben$cross_eff
#' result$M2_ben$e
#' 
#' # Example 3.
#' # Benevolent formulation (method III). Input-oriented.
#' data("Golany_Roll_1989")
#' data_example <- read_data(datadea = Golany_Roll_1989, 
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- cross_efficiency(data_example, 
#'                            orientation = "io", 
#'                            selfapp = TRUE)
#' result$M3_ben$cross_eff
#' result$M3_ben$e
#'   
#' # Example 4.
#' # Arbitrary formulation. Output-oriented.
#' data("Golany_Roll_1989")
#' data_example <- read_data(datadea = Golany_Roll_1989,
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- cross_efficiency(data_example, 
#'                            orientation = "oo", 
#'                            selfapp = TRUE)
#' result$Arbitrary$cross_eff
#' result$Arbitrary$e
#' 
#' # Example 5.
#' # Arbitrary formulation. Input-oriented model under vrs returns-to-scale.
#' data("Lim_Zhu_2015")
#' data_example <- read_data(Lim_Zhu_2015,
#'                           ni = 1, 
#'                           no = 5)
#' cross <- cross_efficiency(data_example,
#'                           epsilon = 0,
#'                           orientation = "io",
#'                           rts = "vrs",
#'                           selfapp = TRUE,
#'                           M2 = FALSE,
#'                           M3 = FALSE)
#' cross$Arbitrary$e
#' 
#' @seealso \code{\link{model_multiplier}}, \code{\link{cross_efficiency_fuzzy}}
#' @export

cross_efficiency <- function(datadea,
                             dmu_eval = NULL,
                             dmu_ref = NULL,
                             epsilon = 0,
                             orientation = c("io", "oo"),
                             rts = c("crs", "vrs"),
                             selfapp = TRUE,
                             correction = FALSE,
                             M2 = TRUE,
                             M3 = TRUE) {
  
  if (is.dea(datadea)) {
    
    deasol <- datadea
    datadea <- deasol$data
    
    if (!any(grepl("multiplier", names(deasol$DMU[[1]])))) {
      stop("Not a multiplier dea solution!")
    }
    
    dmunames <- datadea$dmunames
    nd <- length(dmunames) # number of dmus
    orientation <- deasol$orientation
    dmu_eval <- deasol$dmu_eval
    nde <- length(dmu_eval)
    dmu_ref <- deasol$dmu_ref
    ndr <- length(dmu_ref)
    if (ndr < 3) {
      stop("More than 2 DMUs are needed in dmu_ref.")
    }
    
    rts <- deasol$rts
    if ((rts != "crs") && (rts != "vrs")) {
      stop("Only constant (crs) or variable (vrs) returns to scale.")
    }
    
    epsilon <- deasol$epsilon
    
  } else if (is.deadata(datadea)) {
    
    orientation <- tolower(orientation)
    orientation <- match.arg(orientation)
    
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
    ndr <- length(dmu_ref)
    if (ndr < 3) {
      stop("More than 2 DMUs are needed in dmu_ref.")
    }
    
    # Checking non-controllable or non-discretionary inputs/outputs
    #if ((!is.null(datadea$nc_inputs)) || (!is.null(datadea$nc_outputs))
    #    || (!is.null(datadea$nd_inputs)) || (!is.null(datadea$nd_outputs))) {
    #  deasol <- model_basic(datadea = datadea, dmu_ref = dmu_ref, orientation = orientation, maxslack = FALSE, compute_multiplier = TRUE)
    #} else {
    deasol <- model_multiplier(datadea = datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, orientation = orientation, rts = rts, epsilon = epsilon, compute_lambda = FALSE)
    #}
    
  } else {
    stop("Input should be a dea or deadata class object!")
  }
  
  input <- datadea$input
  output <- datadea$output
  
  eff <- unlist(lapply(deasol$DMU, function(x) x$efficiency))
  mul_input <- do.call(rbind, lapply(deasol$DMU, function(x) x$multiplier_input))
  mul_output <- do.call(rbind, lapply(deasol$DMU, function(x) x$multiplier_output))
  mul_rts <- unlist(lapply(deasol$DMU, function(x) x$multiplier_rts))
  
  if (nrow(mul_input) < nde) {
    stop("There are unfeasible results.")
  }
  
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  if (rts == "crs") {
    if (orientation == "io") {
      cross_eff <- (mul_output %*% output[, dmu_eval] + matrix(mul_rts, nrow = nde, ncol = nde)) /
        mul_input %*% input[, dmu_eval]
    } else {
      cross_eff <-  mul_output %*% output[, dmu_eval] /
        (mul_input %*% input[, dmu_eval] + matrix(mul_rts, nrow = nde, ncol = nde))
    }
  } else {
    if (orientation == "io") {
      if (correction) {
        cross_eff <- mul_output %*% output[, dmu_eval] /
          (mul_input %*% input[, dmu_eval] - matrix(mul_rts, nrow = nde, ncol = nde))
        diag(cross_eff) <- eff
      } else {
        cross_eff <- (mul_output %*% output[, dmu_eval] + matrix(mul_rts, nrow = nde, ncol = nde)) /
        mul_input %*% input[, dmu_eval]
      }
    } else {
      cross_eff <-  mul_output %*% output[, dmu_eval] /
        (mul_input %*% input[, dmu_eval] + matrix(mul_rts, nrow = nde, ncol = nde))
    }
  }
  colnames(cross_eff) <- dmunames[dmu_eval]
  rownames(cross_eff) <- dmunames[dmu_eval]
  
  A <- rep(0, nde)
  e <- A
  
  if (selfapp) {
    A <- rowMeans(cross_eff)
    e <- colMeans(cross_eff)
  } else {
    for (i in 1:nde) {
      A[i] <- sum(cross_eff[i, -i]) / (nde - 1)
      e[i] <- sum(cross_eff[-i, i]) / (nde - 1)
    }
  }
  names(A) <- dmunames[dmu_eval]
  names(e) <- dmunames[dmu_eval]
  
  if (orientation == "io") {
    maverick <- (eff - e) / pmin(eff, e)
  } else {
    maverick <- ((1 / eff) - e) / pmin(1 / eff, e)
  }
  names(maverick) <- dmunames[dmu_eval]
  
  if (orientation == "io") {
    Arbitrary <- list(multiplier_input = mul_input,
                      multiplier_output = mul_output,
                      multiplier_rts = mul_rts,
                      cross_eff = cross_eff,
                      efficiency = eff,
                      e = e,
                      A = A,
                      maverick = maverick)
    orient <- 1
  } else {
    Arbitrary <- list(multiplier_input = mul_input,
                      multiplier_output = mul_output,
                      multiplier_rts = mul_rts,
                      cross_eff = 1 / cross_eff,
                      efficiency = eff,
                      e = 1 / e,
                      A = 1 / A,
                      maverick = maverick)
    input <- -datadea$output
    output <- -datadea$input
    orient <- -1
  }
  
  ########## Aggressive/benevolent computations ##########
  
  M2_agg <- NULL
  M2_ben <- NULL
  M3_agg <- NULL
  M3_ben <- NULL
  
  if (rts == "crs") {
    
    inputnames <- rownames(input)
    outputnames <- rownames(output)
    ni <- nrow(input) # number of  inputs
    no <- nrow(output) # number of outputs
    inputref <- matrix(input[, dmu_ref], nrow = ni) 
    outputref <- matrix(output[, dmu_ref], nrow = no)
    
    mul_input_agg <- matrix(0, nrow = nde, ncol = ni)
    mul_output_agg <- matrix(0, nrow = nde, ncol = no)
    rownames(mul_input_agg) <- dmunames[dmu_eval]
    rownames(mul_output_agg) <- dmunames[dmu_eval]
    colnames(mul_input_agg) <- inputnames
    colnames(mul_output_agg) <- outputnames
    
    mul_input_ben <- mul_input_agg
    mul_output_ben <- mul_output_agg
    
    if (epsilon > 0) {
      f.con.eps <- diag(ni + no)
      f.dir.eps <- rep(">=", ni + no)
      f.rhs.eps <- rep(epsilon, ni + no)
    } else {
      f.con.eps <- NULL
      f.dir.eps <- NULL
      f.rhs.eps <- NULL
    }
    #f.con.eps <- rbind(f.con.eps, diag(ni + no + 2))
    #f.dir.eps <- c(f.dir.eps, rep("<=", ni + no + 2))
    #f.rhs.eps <- c(f.rhs.eps, rep(1e12, ni + no + 2))
    
    ##### Aggressive/benevolent method II #####
    
    if (M2) {
      
      for (i in 1:nde) {
        
        ii <- dmu_eval[i]
        
        if (ii %in% dmu_ref) {
          inputrefnoi <- matrix(inputref[, -which(dmu_ref == ii)], nrow = ni)
          outputrefnoi <- matrix(outputref[, -which(dmu_ref == ii)], nrow = no)
          ndrnoi <- ndr - 1
        } else {
          inputrefnoi <- inputref
          outputrefnoi <- outputref
          ndrnoi <- ndr
        }
        f.dir <- c(rep("<=", ndrnoi), "=", "=", f.dir.eps)
        f.rhs <- c(rep(0, ndrnoi), orient, orient * eff[i], f.rhs.eps)
        #f.dir <- c(rep("<=", ndr), "=", "=", f.dir.eps)
        #f.rhs <- c(rep(0, ndr), orient, orient * eff[i], f.rhs.eps)
        
        f.obj <- c(-rowSums(inputrefnoi), rowSums(outputrefnoi))
        
        f.con.1 <- cbind(-t(inputrefnoi), t(outputrefnoi))
        #f.con.1 <- cbind(-t(inputref), t(outputref))
        f.con.2 <- c(input[, ii], rep(0, no))
        f.con.3 <- c(rep(0, ni), output[, ii])
        f.con <- rbind(f.con.1, f.con.2, f.con.3, f.con.eps)
        
        res <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
        
        mul_input_agg[i, ] <- res[1 : ni]
        mul_output_agg[i, ] <- res[(ni + 1) : (ni + no)]
        
        res <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
        
        mul_input_ben[i, ] <- res[1 : ni]
        mul_output_ben[i, ] <- res[(ni + 1) : (ni + no)]
        
      }
      
      ## Aggressive ##
      
      if (orientation == "io") {
        cross_eff <- mul_output_agg %*% output[, dmu_eval] /
          mul_input_agg %*% input[, dmu_eval]
      } else {
        cross_eff <- mul_input_agg %*% input[, dmu_eval] /
          mul_output_agg %*% output[, dmu_eval]
      }
      
      colnames(cross_eff) <- dmunames[dmu_eval]
      rownames(cross_eff) <- dmunames[dmu_eval]
      
      if (selfapp) {
        A <- rowMeans(cross_eff)
        e <- colMeans(cross_eff)
      } else {
        for (i in 1:nde) {
          A[i] <- sum(cross_eff[i, -i]) / (nde - 1)
          e[i] <- sum(cross_eff[-i, i]) / (nde - 1)
        }
      }
      names(A) <- dmunames[dmu_eval]
      names(e) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        maverick <- (eff - e) / pmin(eff, e)
      } else {
        maverick <- ((1 / eff) - e) / pmin(1 / eff, e)
      }
      names(maverick) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        M2_agg <- list(multiplier_input = mul_input_agg,
                       multiplier_output = mul_output_agg,
                       cross_eff = cross_eff,
                       e = e,
                       A = A,
                       maverick = maverick)
      } else {
        M2_agg <- list(multiplier_input = mul_output_agg,
                       multiplier_output = mul_input_agg,
                       cross_eff = 1 / cross_eff,
                       e =  1 / e,
                       A = 1 / A,
                       maverick = maverick)
      }
      
      ## Benevolent ##
      
      if (orientation == "io") {
        cross_eff <- mul_output_ben %*% output[, dmu_eval] /
          mul_input_ben %*% input[, dmu_eval]
      } else {
        cross_eff <- mul_input_ben %*% input[, dmu_eval] /
          mul_output_ben %*% output[, dmu_eval]
      }
      colnames(cross_eff) <- dmunames[dmu_eval]
      rownames(cross_eff) <- dmunames[dmu_eval]
      
      if (selfapp) {
        A <- rowMeans(cross_eff)
        e <- colMeans(cross_eff)
      } else {
        for (i in 1:nde) {
          A[i] <- sum(cross_eff[i, -i]) / (nde - 1)
          e[i] <- sum(cross_eff[-i, i]) / (nde - 1)
        }
      }
      names(A) <- dmunames[dmu_eval]
      names(e) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        maverick <- (eff - e) / pmin(eff, e)
      } else {
        maverick <- ((1 / eff) - e) / pmin(1 / eff, e)
      }
      names(maverick) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        M2_ben <- list(multiplier_input = mul_input_ben,
                       multiplier_output = mul_output_ben,
                       cross_eff = cross_eff,
                       e = e,
                       A = A,
                       maverick = maverick)
      } else {
        M2_ben <- list(multiplier_input = mul_output_ben,
                       multiplier_output = mul_input_ben,
                       cross_eff = 1 / cross_eff,
                       e = 1 / e,
                       A = 1 / A,
                       maverick = maverick)
      }
      
    }
    
    ##### Aggressive/benevolent method III #####
    
    if (M3) {
      
      for (i in 1:nde) {
        
        ii <- dmu_eval[i]
        
        if (i %in% dmu_ref) {
          inputrefnoi <- matrix(inputref[, -which(dmu_ref == ii)], nrow = ni)
          outputrefnoi <- matrix(outputref[, -which(dmu_ref == ii)], nrow = no)
          ndrnoi <- ndr - 1
        } else {
          inputrefnoi <- inputref
          outputrefnoi <- outputref
          ndrnoi <- ndr
        }
        f.dir <- c(rep("<=", ndrnoi), "=", "=", f.dir.eps)
        f.rhs <- c(rep(0, ndrnoi), orient, 0, f.rhs.eps)
        
        f.obj <- c(rep(0, ni), rowSums(outputrefnoi))
        
        f.con.1 <- cbind(-t(inputrefnoi), t(outputrefnoi))
        f.con.2 <- c(rowSums(inputrefnoi), rep(0, no))
        f.con.3 <- c(-eff[i] * input[, ii], output[, ii])
        f.con <- rbind(f.con.1, f.con.2, f.con.3, f.con.eps)
        
        res <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
        
        mul_input_agg[i, ] <- res[1 : ni]
        mul_output_agg[i, ] <- res[(ni + 1) : (ni + no)]
        
        res <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
        
        mul_input_ben[i, ] <- res[1 : ni]
        mul_output_ben[i, ] <- res[(ni + 1) : (ni + no)]
        
      }
      
      ## Aggressive ##
      
      if (orientation == "io") {
        cross_eff <- mul_output_agg %*% output[, dmu_eval] /
          mul_input_agg %*% input[, dmu_eval]
      } else {
        cross_eff <- mul_input_agg %*% input[, dmu_eval] /
          mul_output_agg %*% output[, dmu_eval]
      }
      colnames(cross_eff) <- dmunames[dmu_eval]
      rownames(cross_eff) <- dmunames[dmu_eval]
      
      if (selfapp) {
        A <- rowMeans(cross_eff)
        e <- colMeans(cross_eff)
      } else {
        for (i in 1:nde) {
          A[i] <- sum(cross_eff[i, -i]) / (nde - 1)
          e[i] <- sum(cross_eff[-i, i]) / (nde - 1)
        }
      }
      names(A) <- dmunames[dmu_eval]
      names(e) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        maverick <- (eff - e) / pmin(eff, e)
      } else {
        maverick <- ((1 / eff) - e) / pmin(1 / eff, e)
      }
      names(maverick) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        M3_agg <- list(multiplier_input = mul_input_agg,
                       multiplier_output = mul_output_agg,
                       cross_eff = cross_eff,
                       e = e,
                       A = A,
                       maverick = maverick)
      } else {
        M3_agg <- list(multiplier_input = mul_output_agg,
                       multiplier_output = mul_input_agg,
                       cross_eff = 1 / cross_eff,
                       e = 1 / e,
                       A = 1 / A,
                       maverick = maverick)
      }
      
      ## Benevolent ##
      
      if (orientation == "io") {
        cross_eff <- mul_output_ben %*% output[, dmu_eval] /
          mul_input_ben %*% input[, dmu_eval]
      } else {
        cross_eff <- mul_input_ben %*% input[, dmu_eval] /
          mul_output_ben %*% output[, dmu_eval]
      }
      colnames(cross_eff) <- dmunames[dmu_eval]
      rownames(cross_eff) <- dmunames[dmu_eval]
      
      if (selfapp) {
        A <- rowMeans(cross_eff)
        e <- colMeans(cross_eff)
      } else {
        for (i in 1:nde) {
          A[i] <- sum(cross_eff[i, -i]) / (nde - 1)
          e[i] <- sum(cross_eff[-i, i]) / (nde - 1)
        }
      }
      names(A) <- dmunames[dmu_eval]
      names(e) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        maverick <- (eff - e) / pmin(eff, e)
      } else {
        maverick <- ((1 / eff) - e) / pmin(1 / eff, e)
      }
      names(maverick) <- dmunames[dmu_eval]
      
      if (orientation == "io") {
        M3_ben <- list(multiplier_input = mul_input_ben,
                       multiplier_output = mul_output_ben,
                       cross_eff = cross_eff,
                       e = e,
                       A = A,
                       maverick = maverick)
      } else {
        M3_ben <- list(multiplier_input = mul_output_ben,
                       multiplier_output = mul_input_ben,
                       cross_eff = 1/ cross_eff,
                       e = 1 / e,
                       A = 1 / A,
                       maverick = maverick)
      }
      
    }
    
  }
  
  ########## Return ##########
  
  deaOutput <- list(orientation = orientation,
              rts = rts,
              selfapp = selfapp,
              Arbitrary = Arbitrary,
              M2_agg = M2_agg,
              M2_ben = M2_ben,
              M3_agg = M3_agg,
              M3_ben = M3_ben,
              data = datadea,
              dmu_eval = dmu_eval,
              dmu_ref = dmu_ref,
              epsilon = epsilon,
              modelname = "cross_efficiency")
  return(structure(deaOutput, class = "dea"))
  
  
}