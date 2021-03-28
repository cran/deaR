#' @title Cross efficiency fuzzy tables
#'   
#' @description Computes the cross-efficiency fuzzy table from dea data or a Guo-Tanaka dea model solution.
#' The (crisp) relative efficiencies for the case \code{h} = 1 are obtained from the CCR model (\code{model_multiplier}).
#' 
#' @usage cross_efficiency_fuzzy(datadea,
#'                        orientation = c("io", "oo"),
#'                        h = 1,
#'                        selfapp = TRUE)
#' 
#' @param datadea An object of class \code{dea_fuzzy} or \code{deadata_fuzzy}. If it is of class \code{dea_fuzzy} it must have
#'                been obtained with \code{modelfuzzy_guotanaka}.
#' @param orientation A string, equal to "io" (input-oriented) or "oo" (output-oriented).
#' @param h A numeric vector with the h-levels (in [0,1]).
#' @param selfapp Logical. If it is \code{TRUE}, self-appraisal is included in the average scores of
#'                \code{A} and \code{e}.
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
#' Doyle, J.; Green, R. (1994). “Efficiency and Cross Efficiency in DEA: Derivations, Meanings and the Uses”,  Journal of Operational Research Society, 45(5), 567–578. \doi{10.2307/2584392}
#'  
#' Guo, P.; Tanaka, H. (2001). "Fuzzy DEA: A Perceptual Evaluation Method", Fuzzy Sets and Systems, 119, 149–160. \doi{10.1016/S0165-0114(99)00106-2}
#' 
#' León, T.; Liern, V.; Ruiz, J.L.; Sirvent, I. (2003). "A Fuzzy Mathematical Programming Approach to the assessment of efficiency with DEA Models", Fuzzy Sets Systems, 139(2), 407–419. \doi{10.1016/S0165-0114(02)00608-5}
#' 
#' Sexton, T.R., Silkman, R.H.; Hogan, A.J. (1986). Data envelopment analysis: critique and extensions. In: Silkman RH (ed) Measuring efficiency: an assessment of data envelopment analysis, vol 32. Jossey-Bass, San Francisco, pp 73–104. \doi{10.1002/ev.1441}  
#' 
#' @examples 
#'  data("Guo_Tanaka_2001")
#'  datadea <- read_data_fuzzy(datadea = Guo_Tanaka_2001, 
#'                             inputs.mL = 2:3, 
#'                             inputs.dL = 4:5, 
#'                             outputs.mL = 6:7, 
#'                             outputs.dL = 8:9)
#'  result <- cross_efficiency_fuzzy(datadea = datadea, 
#'                                   h = seq(0, 1, 0.2))
#'  
#' @export

cross_efficiency_fuzzy <- function(datadea,
                                   orientation = c("io", "oo"),
                                   h = 1,
                                   selfapp = TRUE) {
  
  if (is.dea_fuzzy(datadea)) {
    
    deasol <- datadea
    if (deasol$modelname != "fuzzy_guotanaka") {
      stop("Not a modelfuzzy_guotanaka dea solution!")
    }
    datadea <- deasol$data
    
    dmunames <- datadea$dmunames
    nd <- length(dmunames) # number of dmus
    if (length(deasol$dmu_eval) != nd) {
      stop("We need multipliers for all DMUs.")
    }
    
    orientation <- deasol$orientation
    h <- deasol$h
    
  } else if (is.deadata_fuzzy(datadea)) {
    
    orientation <- tolower(orientation)
    orientation <- match.arg(orientation)
    
    dmunames <- datadea$dmunames
    nd <- length(dmunames) # number of dmus
    
    deasol <- modelfuzzy_guotanaka(datadea = datadea, orientation = orientation, h = h)
    
  } else {
    stop("Input should be a dea_fuzzy or deadata_fuzzy class object!")
  }
  
  nh <- length(h) # number of h-levels
  hlevel <- vector(mode = "list", length = nh)
  names(hlevel) <- as.character(h)
  
  if (orientation == "io") {
    input.d <- datadea$input$dL # spread
    orient <- 1
    objg0 <- "max"
    g0.dir <- c("=", "<=")
    f.dir <- c("=", ">=", g0.dir, rep("<=", nd + nd))
    #f.dir <- c("=", ">=", g0.dir, rep("<=", nd + nd - 2))
    #f.dir <- c("=", ">=", g0.dir, rep("<=", nd + nd), rep("=", nd - 1))
  } else {
    input.d <- -datadea$output$dL # spread
    orient <- -1
    objg0 <- "min"
    g0.dir <- c("=", ">=")
    f.dir <- c("=", "<=", g0.dir, rep("<=", nd + nd))
    #f.dir <- c("=", "<=", g0.dir, rep("<=", nd + nd - 2))
    #f.dir <- c("=", "<=", g0.dir, rep("<=", nd + nd), rep("=", nd - 1))
  }
  
  for (i in 1:nh) {
    
    input.m <- datadea$input$mL # center
    output.m <- datadea$output$mL
    inputnames <- rownames(input.m)
    outputnames <- rownames(output.m)
    ni <- nrow(input.m) # number of  inputs
    no <- nrow(output.m) # number of outputs
    
    input.L <- deasol$hlevel[[i]]$input$Lower
    input.U <- deasol$hlevel[[i]]$input$Upper
    output.L <- deasol$hlevel[[i]]$output$Lower
    output.U <- deasol$hlevel[[i]]$output$Upper
    mul_input <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x) x$multiplier_input))
    mul_output <- do.call(rbind, lapply(deasol$hlevel[[i]]$DMU, function(x) x$multiplier_output))
    
    cross_eff.m <- (mul_output %*% output.m) / (mul_input %*% input.m)
    cross_eff.dL <- cross_eff.m - (mul_output %*% output.L) / (mul_input %*% input.U)
    cross_eff.dR <- (mul_output %*% output.U) / (mul_input %*% input.L) - cross_eff.m
    
    colnames(cross_eff.m) <- dmunames
    rownames(cross_eff.m) <- dmunames
    colnames(cross_eff.dL) <- dmunames
    rownames(cross_eff.dL) <- dmunames
    colnames(cross_eff.dR) <- dmunames
    rownames(cross_eff.dR) <- dmunames
    
    eff.m <- diag(cross_eff.m)
    eff.dL <- diag(cross_eff.dL)
    eff.dR <- diag(cross_eff.dR)
    
    A.m <- rep(0, nd)
    names(A.m) <- dmunames
    A.dL <- A.m
    A.dR <- A.m
    e.m <- A.m
    e.dL <- A.m
    e.dR <- A.m
    
    if (selfapp) {
      A.m <- rowMeans(cross_eff.m)
      A.dL <- rowMeans(cross_eff.dL)
      A.dR <- rowMeans(cross_eff.dR)
      e.m <- colMeans(cross_eff.m)
      e.dL <- colMeans(cross_eff.dL)
      e.dR <- colMeans(cross_eff.dR)
    } else {
      for (j in 1:nd) {
        A.m[j] <- sum(cross_eff.m[j, -j]) / (nd - 1)
        A.dL[j] <- sum(cross_eff.dL[j, -j]) / (nd - 1)
        A.dR[j] <- sum(cross_eff.dR[j, -j]) / (nd - 1)
        e.m[j] <- sum(cross_eff.m[-j, j]) / (nd - 1)
        e.dL[j] <- sum(cross_eff.dL[-j, j]) / (nd - 1)
        e.dR[j] <- sum(cross_eff.dR[-j, j]) / (nd - 1)
      }
    }
    
    #maverick.m <- (eff.m - e.m) / pmin(eff.m, e.m)
    #maverick.dL <- maverick.m - (eff.m - eff.dL - e.m  - e.dR) / pmin(eff.m + eff.dR, e.m + e.dR)
    #maverick.dR <- (eff.m + eff.dR - e.m  + e.dL) / pmin(eff.m - eff.dL, e.m - e.dL) - maverick.m
    #names(maverick.m) <- dmunames
    #names(maverick.dL) <- dmunames
    #names(maverick.dR) <- dmunames
    
    Y2 = e.m + 0.25 * (e.dR - e.dL)
    names(Y2) <- dmunames
    
    if (orientation == "io") {
      Arbitrary <- list(multiplier_input = mul_input,
                        multiplier_output = mul_output,
                        cross_eff = list(dL = cross_eff.dL, m = cross_eff.m, dR = cross_eff.dR),
                        e = list(dL = e.dL, m = e.m, dR = e.dR),
                        A = list(dL = A.dL, m = A.m, dR = A.dR),
                        #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                        Y2 = Y2)
    } else {
      Arbitrary <- list(multiplier_input = mul_input,
                        multiplier_output = mul_output,
                        cross_eff = list(dL = (1 / cross_eff.m) - (1 / (cross_eff.m + cross_eff.dR)), m = 1 / cross_eff.m, dR =  (1 / (cross_eff.m - cross_eff.dL)) - (1 / cross_eff.m)),
                        e = list(dL = (1 / e.m) - (1 / (e.m + e.dR)), m = 1 / e.m, dR =  (1 / (e.m - e.dL)) - (1 / e.m)),
                        A = list(dL = (1 / A.m) - (1 / (A.m + A.dR)), m = 1 / A.m, dR =  (1 / (A.m - A.dL)) - (1 / A.m)),
                        #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                        Y2 = Y2)
      
    }
    
    ########## Aggressive/benevolent computations ##########
    
    if (orientation == "oo") {
      input.m <- -datadea$output$mL
      output.m <- -datadea$input$mL
      input.L <- -deasol$hlevel[[i]]$output$Lower
      input.U <- -deasol$hlevel[[i]]$output$Upper
      output.L <- -deasol$hlevel[[i]]$input$Lower
      output.U <- -deasol$hlevel[[i]]$input$Upper
      mul_output <- mul_input
    }
    inputnames <- rownames(input.L)
    outputnames <- rownames(output.L)
    ni <- nrow(input.L) # number of  inputs
    no <- nrow(output.L) # number of outputs
    
    a <- h[i]
    e_spread <- max(input.d / input.m)
    one.L <- 1 - (1 - a) * e_spread
    one.U <- 1 + (1 - a) * e_spread
    
    mul_input_agg <- matrix(0, nrow = nd, ncol = ni)
    mul_output_agg <- matrix(0, nrow = nd, ncol = no)
    rownames(mul_input_agg) <- dmunames
    rownames(mul_output_agg) <- dmunames
    colnames(mul_input_agg) <- inputnames
    colnames(mul_output_agg) <- outputnames
    
    mul_input_ben <- mul_input_agg
    mul_output_ben <- mul_output_agg
    
    g0.rhs <- orient * c(one.L, one.U)
    
    for (j in 1:nd) {
      
      beta <- mul_output[j, ] %*% output.L[, j] # objval of original program
      
      # Compute g0
      g0.obj <- input.d[, j]
      g0.con <- rbind(input.L[, j], input.U[, j])
      g0 <- lp(objg0, g0.obj, g0.con, g0.dir, g0.rhs)$objval
      
      # Vector de coeficientes de la función objetivo
      f.obj <- c(rowSums(matrix(input.U[, -j], nrow = ni)), -rowSums(matrix(output.L[, -j], nrow = no))) 
      #f.obj <- c(rep(0, ni + no), rep(1, nd - 1))
      
      # Matriz técnica
      f.con.0 <- c(rep(0, ni), output.L[, j])
      f.con.1 <- c(g0.obj, rep(0, no))
      f.con.2 <- cbind(g0.con, matrix(0, nrow = 2, ncol = no))
      f.con.3 <- cbind(-t(input.L), t(output.L))
      f.con.4 <- cbind(-t(input.U), t(output.U))
      #f.con.3 <- cbind(-t(input.L[, -j]), t(output.L[, -j]))
      #f.con.4 <- cbind(-t(input.U[, -j]), t(output.U[, -j]))
      f.con <- rbind(f.con.0, f.con.1, f.con.2, f.con.3, f.con.4)
      #f.con.0 <- c(rep(0, ni), output.L[, j], rep(0, nd - 1))
      #f.con.1 <- c(g0.obj, rep(0, no), rep(0, nd - 1))
      #f.con.2 <- cbind(g0.con, matrix(0, nrow = 2, ncol = no), matrix(0, nrow = 2, ncol = nd - 1))
      #f.con.3 <- cbind(-t(input.L), t(output.L), matrix(0, nrow = nd, ncol = nd - 1))
      #f.con.4 <- cbind(-t(input.U), t(output.U), matrix(0, nrow = nd, ncol = nd - 1))
      #f.con.5 <- cbind(-t(input.U[, -j]), t(output.L[, -j]), orient * diag(nd - 1))
      #f.con <- rbind(f.con.0, f.con.1, f.con.2, f.con.3, f.con.4, f.con.5)
      
      # Vector de términos independientes
      f.rhs <- c(beta, g0, orient * c(one.L, one.U), rep(0, nd + nd))
      #f.rhs <- c(beta, g0, orient * c(one.L, one.U), rep(0, nd + nd - 2))
      #f.rhs <- c(beta, g0, orient * c(one.L, one.U), rep(0, nd + nd + nd - 1))
      
      res <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
      
      mul_input_agg[j, ] <- res[1 : ni]
      mul_output_agg[j, ] <- res[(ni + 1) : (ni + no)]
      
      res <- lp("min", f.obj, f.con, f.dir, f.rhs)$solution
      
      mul_input_ben[j, ] <- res[1 : ni]
      mul_output_ben[j, ] <- res[(ni + 1) : (ni + no)]
      
    }
    
    ## Aggressive ##
    
    if (orientation == "io") {
      cross_eff.m <- (mul_output_agg %*% output.m) / (mul_input_agg %*% input.m)
      cross_eff.dL <- cross_eff.m - (mul_output_agg %*% output.L) / (mul_input_agg %*% input.U)
      cross_eff.dR <- (mul_output_agg %*% output.U) / (mul_input_agg %*% input.L) - cross_eff.m
    } else {
      cross_eff.m <- (mul_input_agg %*% input.m) / (mul_output_agg %*% output.m)
      cross_eff.dL <- cross_eff.m - (mul_input_agg %*% input.L) / (mul_output_agg %*% output.U)
      cross_eff.dR <- (mul_input_agg %*% input.U) / (mul_output_agg %*% output.L) - cross_eff.m
    }
    
    colnames(cross_eff.m) <- dmunames
    rownames(cross_eff.m) <- dmunames
    colnames(cross_eff.dL) <- dmunames
    rownames(cross_eff.dL) <- dmunames
    colnames(cross_eff.dR) <- dmunames
    rownames(cross_eff.dR) <- dmunames
    
    eff.m <- diag(cross_eff.m)
    eff.dL <- diag(cross_eff.dL)
    eff.dR <- diag(cross_eff.dR)
    
    A.m <- rep(0, nd)
    names(A.m) <- dmunames
    A.dL <- A.m
    A.dR <- A.m
    e.m <- A.m
    e.dL <- A.m
    e.dR <- A.m
    
    if (selfapp) {
      A.m <- rowMeans(cross_eff.m)
      A.dL <- rowMeans(cross_eff.dL)
      A.dR <- rowMeans(cross_eff.dR)
      e.m <- colMeans(cross_eff.m)
      e.dL <- colMeans(cross_eff.dL)
      e.dR <- colMeans(cross_eff.dR)
    } else {
      for (j in 1:nd) {
        A.m[j] <- sum(cross_eff.m[j, -j]) / (nd - 1)
        A.dL[j] <- sum(cross_eff.dL[j, -j]) / (nd - 1)
        A.dR[j] <- sum(cross_eff.dR[j, -j]) / (nd - 1)
        e.m[j] <- sum(cross_eff.m[-j, j]) / (nd - 1)
        e.dL[j] <- sum(cross_eff.dL[-j, j]) / (nd - 1)
        e.dR[j] <- sum(cross_eff.dR[-j, j]) / (nd - 1)
      }
    }
    
    #maverick.m <- (eff.m - e.m) / pmin(eff.m, e.m)
    #maverick.dL <- maverick.m - (eff.m - eff.dL - e.m  - e.dR) / pmin(eff.m + eff.dR, e.m + e.dR)
    #maverick.dR <- (eff.m + eff.dR - e.m  + e.dL) / pmin(eff.m - eff.dL, e.m - e.dL) - maverick.m
    #names(maverick.m) <- dmunames
    #names(maverick.dL) <- dmunames
    #names(maverick.dR) <- dmunames
    
    Y2 = e.m + 0.25 * (e.dR - e.dL)
    names(Y2) <- dmunames
    
    if (orientation == "io") {
      Agg <- list(multiplier_input = mul_input_agg,
                  multiplier_output = mul_output_agg,
                  cross_eff = list(dL = cross_eff.dL, m = cross_eff.m, dR = cross_eff.dR),
                  e = list(dL = e.dL, m = e.m, dR = e.dR),
                  A = list(dL = A.dL, m = A.m, dR = A.dR),
                  #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                  Y2 = Y2)
    } else {
      Agg <- list(multiplier_input = mul_output_agg,
                  multiplier_output = mul_input_agg,
                  cross_eff = list(dL = (1 / cross_eff.m) - (1 / (cross_eff.m + cross_eff.dR)), m = 1 / cross_eff.m, dR =  (1 / (cross_eff.m - cross_eff.dL)) - (1 / cross_eff.m)),
                  e = list(dL = (1 / e.m) - (1 / (e.m + e.dR)), m = 1 / e.m, dR =  (1 / (e.m - e.dL)) - (1 / e.m)),
                  A = list(dL = (1 / A.m) - (1 / (A.m + A.dR)), m = 1 / A.m, dR =  (1 / (A.m - A.dL)) - (1 / A.m)),                  
                  #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                  Y2 = Y2)
    }
    
    ## Benevolent ##
    
    if (orientation == "io") {
      cross_eff.m <- (mul_output_ben %*% output.m) / (mul_input_ben %*% input.m)
      cross_eff.dL <- cross_eff.m - (mul_output_ben %*% output.L) / (mul_input_ben %*% input.U)
      cross_eff.dR <- (mul_output_ben %*% output.U) / (mul_input_ben %*% input.L) - cross_eff.m
    } else {
      cross_eff.m <- (mul_input_ben %*% input.m) / (mul_output_ben %*% output.m)
      cross_eff.dL <- cross_eff.m - (mul_input_ben %*% input.L) / (mul_output_ben %*% output.U)
      cross_eff.dR <- (mul_input_ben %*% input.U) / (mul_output_ben %*% output.L) - cross_eff.m
    }
    
    colnames(cross_eff.m) <- dmunames
    rownames(cross_eff.m) <- dmunames
    colnames(cross_eff.dL) <- dmunames
    rownames(cross_eff.dL) <- dmunames
    colnames(cross_eff.dR) <- dmunames
    rownames(cross_eff.dR) <- dmunames
    
    eff.m <- diag(cross_eff.m)
    eff.dL <- diag(cross_eff.dL)
    eff.dR <- diag(cross_eff.dR)
    
    A.m <- rep(0, nd)
    names(A.m) <- dmunames
    A.dL <- A.m
    A.dR <- A.m
    e.m <- A.m
    e.dL <- A.m
    e.dR <- A.m
    
    if (selfapp) {
      A.m <- rowMeans(cross_eff.m)
      A.dL <- rowMeans(cross_eff.dL)
      A.dR <- rowMeans(cross_eff.dR)
      e.m <- colMeans(cross_eff.m)
      e.dL <- colMeans(cross_eff.dL)
      e.dR <- colMeans(cross_eff.dR)
    } else {
      for (j in 1:nd) {
        A.m[j] <- sum(cross_eff.m[j, -j]) / (nd - 1)
        A.dL[j] <- sum(cross_eff.dL[j, -j]) / (nd - 1)
        A.dR[j] <- sum(cross_eff.dR[j, -j]) / (nd - 1)
        e.m[j] <- sum(cross_eff.m[-j, j]) / (nd - 1)
        e.dL[j] <- sum(cross_eff.dL[-j, j]) / (nd - 1)
        e.dR[j] <- sum(cross_eff.dR[-j, j]) / (nd - 1)
      }
    }
    
    #maverick.m <- (eff.m - e.m) / pmin(eff.m, e.m)
    #maverick.dL <- maverick.m - (eff.m - eff.dL - e.m  - e.dR) / pmin(eff.m + eff.dR, e.m + e.dR)
    #maverick.dR <- (eff.m + eff.dR - e.m  + e.dL) / pmin(eff.m - eff.dL, e.m - e.dL) - maverick.m
    #names(maverick.m) <- dmunames
    #names(maverick.dL) <- dmunames
    #names(maverick.dR) <- dmunames
    
    Y2 = e.m + 0.25 * (e.dR - e.dL)
    names(Y2) <- dmunames
    
    if (orientation == "io") {
      Ben <- list(multiplier_input = mul_input_ben,
                  multiplier_output = mul_output_ben,
                  cross_eff = list(dL = cross_eff.dL, m = cross_eff.m, dR = cross_eff.dR),
                  e = list(dL = e.dL, m = e.m, dR = e.dR),
                  A = list(dL = A.dL, m = A.m, dR = A.dR),
                  #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                  Y2 = Y2)
    } else {
      Ben <- list(multiplier_input = mul_output_ben,
                  multiplier_output = mul_input_ben,
                  cross_eff = list(dL = (1 / cross_eff.m) - (1 / (cross_eff.m + cross_eff.dR)), m = 1 / cross_eff.m, dR =  (1 / (cross_eff.m - cross_eff.dL)) - (1 / cross_eff.m)),
                  e = list(dL = (1 / e.m) - (1 / (e.m + e.dR)), m = 1 / e.m, dR =  (1 / (e.m - e.dL)) - (1 / e.m)),
                  A = list(dL = (1 / A.m) - (1 / (A.m + A.dR)), m = 1 / A.m, dR =  (1 / (A.m - A.dL)) - (1 / A.m)),                  
                  #maverick = list(dL = maverick.dL, m = maverick.m, dR = maverick.dR),
                  Y2 = Y2)
    }
    
    hlevel[[i]] <- list(input = deasol$hlevel[[i]]$input,
                        output = deasol$hlevel[[i]]$output,
                        Arbitrary = Arbitrary,
                        Agg = Agg,
                        Ben = Ben)
    
  }
  
  ########## Return ##########
  
  return(list(orientation = orientation,
              rts = "crs",
              selfapp = selfapp,
              h = h,
              hlevel = hlevel,
              data = datadea))
  
}