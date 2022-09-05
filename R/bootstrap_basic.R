#' @title Bootstrapping DEA
#'   
#' @description To bootstrap efficiency scores, deaR uses the algorithm proposed by Simar and Wilson (1998). For now, the function bootstrap_basic can only be used with basic DEA models.
#' 
#' @usage bootstrap_basic(datadea,
#'                 orientation = c("io", "oo"),
#'                 rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'                 L = 1,
#'                 U = 1,
#'                 B = 2000,
#'                 h = NULL,
#'                 alpha = 0.05)
#' 
#' @param datadea The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param orientation A string, equal to "io" (input oriented) or "oo" (output oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#'            "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param B Number of bootstrap iterations.
#' @param h Bandwidth of smoothing window. By default h=0.014 (You can set h to any other value). The optimal bandwidth factor can also be calculated following the proposals of Silverman (1986) and Daraio y Simar (2007). So, h="h1" is the optimal h referred as "robust normal-reference rule" (Daraio and Simar, 2007 p.60), h="h2" is the value of h1 but instead of the factor 1.06 with the factor 0.9, h="h3" is the value of h1 adjusted for scale and sample size (Daraio and Simar, 2007 p.61), and h="h4" is the bandwidth provided by a Gaussian kernel density estimate. 
#' @param alpha Between 0 and 1 (for confidence intervals).
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
#' Behr, A. (2015). Production and Efficiency Analysis with R. Springer.
#' 
#' Bogetoft, P.; Otto, L. (2010).  Benchmarking with DEA, SFA, and R. Springer.
#' 
#' Daraio, C.; Simar, L. (2007). Advanced Robust and Nonparametric Methods in Efficiency Analysis: Methodology and Applications. New York: Springer.
#' 
#' Färe, R.; Grosskopf, S.; Kokkenlenberg, E. (1989). "Measuring Plant Capacity, Utilization and Technical Change: A Nonparametric Approach". International Economic Review, 30(3), 655-666. 
#' 
#' Löthgren, M.; Tambour, M. (1999). "Bootstrapping the Data Envelopment Analysis Malmquist Productivity Index". Applied Economics, 31, 417-425. 
#' 
#' Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis. London: Chapman and Hall. 
#' 
#' Simar, L.; Wilson, P.W. (1998). "Sensitivity Analysis of Efficiency Scores: How to Bootstrap in Nonparametric Frontier Models". Management Science, 44(1), 49-61. 
#' 
#' Simar, L.; Wilson, P.W. (1999). "Estimating and Bootstrapping Malmquist Indices". European Journal of Operational Research, 115, 459-471. 
#' 
#' Simar, L.; Wilson, P.W. (2008). Statistical Inference in Nonparametric Frontier Models: Recent Developments and Perspective. In H.O. Fried; C.A. Knox Lovell and S.S. Schmidt (eds.) The Measurement of Productive Efficiency and Productivity Growth. New York: Oxford University Press. \doi{10.1093/acprof:oso/9780195183528.001.0001}
#' 
#' @examples 
#' # To replicate the results in Simar y Wilson (1998, p. 58) you have to
#' # set B=2000 (in the example B = 100 to save time)
#' data("Electric_plants")
#' data_example <- read_data(Electric_plants, 
#'                           ni = 3, 
#'                           no = 1)
#' result <- bootstrap_basic(datadea = data_example,
#'                              orientation = "io",
#'                              rts = "vrs",
#'                              B = 100)
#' result$score_bc
#' result$CI
#' 
#' @import lpSolve 
#' @importFrom stats IQR median quantile rnorm runif sd var density
#' 
#' @export
  
bootstrap_basic <- function(datadea,
                            orientation = c("io", "oo"),
                            rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
                            L = 1,
                            U = 1,
                            B = 2000,
                            h = NULL,
                            alpha = 0.05) {
 
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
  if (!is.null(datadea$nc_inputs) || !is.null(datadea$nc_outputs)) {
    warning("This model does not take into account the non-controllable feature for inputs/outputs.")
  }
  if (!is.null(datadea$nd_inputs) || !is.null(datadea$nd_outputs)) {
    warning("This model does not take into account the non-discretionary feature for inputs/outputs.")
  }
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
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
  nd <- length(dmunames)
  if (orientation == "io") {
    input <- datadea$input
    output <- datadea$output
    obj <- "min"
  } else {
    input <- -datadea$output
    output <- -datadea$input
    obj <- "max"
  }
  
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input)
  no <- nrow(output)
  score <- rep(0, nd)
  names(score) <- dmunames
  score_bc <- score
  f.obj <- c(1, rep(0, nd))
  if (rts == "crs") {
    f.con.rs <- NULL
    f.con2.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- cbind(0, matrix(1, nrow = 1, ncol = nd))
    f.con2.rs <- cbind(matrix(1, nrow = 1, ncol = nd), matrix(0, 
                                                              nrow = 1, ncol = ni + no))
    f.rhs.rs <- 1
    if (rts == "vrs") {
      f.dir.rs <- "="
    }     else if (rts == "nirs") {
      f.dir.rs <- "<="
    }    else if (rts == "ndrs") {
      f.dir.rs <- ">="
    }    else {
      f.con.rs <- rbind(f.con.rs, f.con.rs)
      f.con2.rs <- rbind(f.con2.rs, f.con2.rs)
      f.dir.rs <- c(">=", "<=")
      f.rhs.rs <- c(L, U)
    }
  }
  f.dir <- c(rep("<=", ni), rep(">=", no), f.dir.rs)
  f.con.2 <- cbind(matrix(0, nrow = no, ncol = 1), output)
  for (i in 1:nd) {
    f.con.1 <- cbind(-input[, i], input)
    f.con <- rbind(f.con.1, f.con.2, f.con.rs)
    f.rhs <- c(rep(0, ni), output[, i], f.rhs.rs)
    score[i] <- lp(obj, f.obj, f.con, f.dir, f.rhs)$solution[1]
  }
  if (orientation == "io") {
    score_sp <- 1/score
  } else {
    score_sp <-  score
  }
  
  new_set <- round(c(score_sp, 2 - score_sp),6)
  
  if (is.null(h)) {
    
    h <- 0.014
    
  }   else if (!is.numeric(h)) {
    if (h %in% c("h1", "h2", "h3", "h4")) {
      sd_new_set <- sd(new_set)
      iqr_new_set <- IQR(new_set)
      if (sd_new_set > iqr_new_set/1.34) {
        desviation <- iqr_new_set/1.34
      }       else {
        desviation <- sd_new_set
      }
      
      if (h == "h1") {
        h <- 1.06 * desviation * length(new_set)^(-1/5)
      }
      else if (h == "h2") {
        h <- 0.9 * desviation * length(new_set)^(-1/5)
      }
      else if (h == "h3") {
        #ajustar hopt (eq. 3.26 de Daraio y Simar (2007) y Simar y Wilson 2006a como en Bogetoft and Otto (2010)
        h <- (0.9 * desviation * length(new_set)^(-1/5)) * 
          (length(new_set)/nd)^(1/5) * (sd(score_sp)/sd(new_set))
      }
      else {
        # bandwidth provided by a Gaussian kernel density estimate
        dens <-  density(new_set)
        h <- dens$bw 
      }
      
    }
    else {
      stop("Incorrect bandwidth argument h.")
    }
  }
  
  estimates_bootstrap <- matrix(nrow = B, ncol = nd)
  colnames(estimates_bootstrap) <- dmunames
  sigma_score <- var(new_set)

  for (b in 1:B) {
    beta <- sample(new_set, nd, replace = TRUE)
    epsilon <- rnorm(nd)
    thetatilde <- beta + h * epsilon
    
    thetatilde_est = 1 + 1/sqrt(1+h^2/sigma_score)*(thetatilde-1)
    
    thetatilde_est <- ifelse(thetatilde_est >= 1,
                             thetatilde_est,
                             2 - thetatilde_est)
    
    if(orientation=="io"){
      input_corrected <- t((thetatilde_est/score_sp) * t(input))  #cambia en funcion de la orientacion
    } else{
      input_corrected <- t((score_sp/thetatilde_est) * t(input))  #cambia en funcion de la orientacion
    }
    
    for (i in 1:nd) {
      f.con.1 <- cbind(-input[, i], input_corrected)
      f.con <- rbind(f.con.1, f.con.2, f.con.rs)
      f.rhs <- c(rep(0, ni), output[, i], f.rhs.rs)
      estimates_bootstrap[b, i] <- lp(obj, f.obj, f.con, 
                                      f.dir, f.rhs)$solution[1]
    }
    
  }
  
  if (orientation == "io") {
    estimates_bootstrap <- 1/estimates_bootstrap
  }
  
  mean_estimates_boot <- apply(estimates_bootstrap, 2, mean, na.rm=TRUE)
  var_estimates_boot <- apply(estimates_bootstrap, 2, var, na.rm=TRUE)
  median_estimates_boot <- apply(estimates_bootstrap, 2, median, na.rm=TRUE)
  bias <- mean_estimates_boot - score_sp
  score_bc <- score_sp - bias
  
  # Intervalos de confianza
  estimates_boot_corrected <- t(apply(-estimates_bootstrap, 1, function(x) x + 2*score_sp))
  CI_low <- apply(estimates_boot_corrected, 2, quantile, alpha/2, na.rm=TRUE)
  CI_up <- apply(estimates_boot_corrected, 2, quantile, 1 - alpha/2, na.rm=TRUE)
  CI <- data.frame(CI_low, CI_up)
  
  # Resultados en input oriented
  if (orientation == "io") {
    score <- 1/score_sp
    estimates_bootstrap <- 1/estimates_bootstrap
    score_bc <- 1/score_bc
    CI <- t(apply(1/CI,1,rev))
    colnames(CI) <- c("CI_low","CI_up")
  }
  
  names(score_bc) <- dmunames
  names(mean_estimates_boot) <- dmunames
  names(var_estimates_boot) <- dmunames
  names(median_estimates_boot) <- dmunames
  names(CI_low) <- dmunames
  names(CI_up) <- dmunames
  descriptives <- data.frame(mean_estimates_boot, var_estimates_boot, 
                             median_estimates_boot)
  
  res <- list(modelname = "bootstrap", orientation = orientation, 
              rts = rts, L = L, U = U, score = score, bandwith = h, score_bc = score_bc, bias = bias,
              descriptives = descriptives, CI = CI, estimates_bootstrap = estimates_bootstrap, 
              data = datadea)
  return(structure(res, class = "dea"))
}
