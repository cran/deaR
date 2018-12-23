#' @title Malmquist index
#'   
#' @description This function calculates the conventional input/output oriented Malmquist index under constant and variables returns-to-scale.
#' 
#' @note In the results: EC=Efficiency Change, PTEC=Pure Technical Efficiency Change, SEC=Scale Efficiency Change, TC= Technological Change, MI= Malmquist Index 
#' @usage malmquist_index(datadealist,
#'                        dmu_eval = NULL,
#'                        dmu_ref = NULL,
#'                        orientation = c("io", "oo"))
#' 
#' @param datadealist A list with the data at different times, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' @param orientation A string, equal to "io" (input oriented) or "oo" (output oriented).
#'   
#' @return A numeric list with Malmquist index and efficiencies.
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
#' Caves, D.W.; Christensen, L. R. y Diewert, W.E. (1982). “The Economic Theory of Index Numbers and the Measurement of Input, Output, and Productivity”. Econometrica, 50(6), 1393-1414. 
#'  
#' Färe, R.; Grosskopf, S.; Lindgren, B.; Roos, P. (1989). “Productivity Developments in Swedish Hospitals: A Malmquist Output Index Approach”. Discussion paper nº 89-3. Southern Illinois University. Illinois.
#' 
#' Färe, R.; Grosskopf, S.; Lindgren, B.; Roos, P. (1992). “Productivity changes in Swedish Pharmacies 1980-89: A nonparametric Malmquist Approach”. Journal of productivity Analysis, 3(3), 85-101. 
#' 
#' Färe, R.; Grosskopf, S.; Norris, M. y Zhang, Z. (1994). “Productivity Growth, Technical Progress, and Efficiency Change in Industrialized Countries”. American Economic Review, 84(1), 66-83. 

#' @examples 
#' # Example 1. With dataset in wide format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("Economy")
#' data_example <- read_malmquist(datadea = Economy,
#'                                nper = 5, 
#'                                arrangement = "horizontal",
#'                                ni = 2, 
#'                                no = 1)
#' result <- malmquist_index(data_example, orientation = "io")
#' mi <- result$mi
#' effch <- result$ec
#' tech <- result$tc
#' 
#' # Example 2. With dataset in long format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("EconomyLong")
#' data_example2 <- read_malmquist(EconomyLong,
#'                                 percol = 2, 
#'                                 arrangement = "vertical",
#'                                 inputs = 3:4, 
#'                                 outputs = 5)
#' result2 <- malmquist_index(data_example2, orientation = "io")
#' mi2 <- result2$mi
#' effch2 <- result2$ec
#' tech2 <- result2$tc
#' 
#' @import lpSolve
#' 
#' @export
  
malmquist_index <- function(datadealist,
                            dmu_eval = NULL,
                            dmu_ref = NULL,
                            orientation = c("io", "oo")) {
  
  nt <- length(datadealist)
  
  if (nt < 2) {
    stop("Number of times should be >= 2.")
  }
  
  # Cheking data
  for (t in 1:nt) {
    if (!is.deadata(datadealist[[t]])) {
      stop("Data should be of class deadata. Run read_data function first!")
    }
  }
  
  dmunames <- datadealist[[1]]$dmunames
  nd <- length(dmunames) # number of dmus
  for (t in 2:nt) {
    if (nd != length(datadealist[[t]]$dmunames)) {
      stop("Number of DMUs does not coincide.")
    }
  }
  
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
  
  ni <- nrow(datadealist[[1]]$input) # number of inputs
  no <- nrow(datadealist[[1]]$output) # number of outputs
  for (t in 2:nt) {
    if (ni != nrow(datadealist[[t]]$input)) {
      stop("Number of inputs does not coincide.")
    }
    if (no != nrow(datadealist[[t]]$output)) {
      stop("Number of outputs does not coincide.")
    }
  }
  
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
  
  mi <- matrix(0, nrow = nt - 1, ncol = nde)
  colnames(mi) <- dmunames[dmu_eval]
  eff <- matrix(0, nrow = nt, ncol = nde)
  colnames(eff) <- dmunames[dmu_eval]
  effvrs <- eff
  eff12 <- mi
  eff21 <- mi
  
  if (orientation == "io") {
    
    obj <- "min"
    f.obj <- c(1, rep(0, ndr))
    
    f.dir <- c(rep("<=", ni), rep(">=", no))
    f.con.vrs <- cbind(0, matrix(1, nrow = 1, ncol = ndr))
    
    for (i in 1:nde) {
      
      ii <- dmu_eval[i]
      
      f.con.1 <- cbind(-datadealist[[1]]$input[, ii], 
                       matrix(datadealist[[1]]$input[, dmu_ref], nrow = ni))
      f.con.2 <- cbind(matrix(0, nrow = no, ncol = 1),
                       matrix(datadealist[[1]]$output[, dmu_ref], nrow = no))
      f.con1 <- rbind(f.con.1, f.con.2)
      
      f.rhs1 <- c(rep(0, ni), datadealist[[1]]$output[, ii])
      
      eff[1, i] <- lp(obj, f.obj, f.con1, f.dir, f.rhs1)$objval
      effvrs[1, i] <- lp(obj, f.obj, rbind(f.con1, f.con.vrs), c(f.dir, "="), c(f.rhs1, 1))$objval
      
      for (t in 2:nt) {
        
        f.con.1 <- cbind(-datadealist[[t]]$input[, ii],
                         matrix(datadealist[[t]]$input[, dmu_ref], nrow = ni))
        f.con.2 <- cbind(matrix(0, nrow = no, ncol = 1),
                         matrix(datadealist[[t]]$output[, dmu_ref], nrow = no))
        f.con2 <- rbind(f.con.1, f.con.2)
        
        f.rhs2 <- c(rep(0, ni), datadealist[[t]]$output[, ii])
        
        eff[t, i] <- lp(obj, f.obj, f.con2, f.dir, f.rhs2)$objval
        effvrs[t, i] <- lp(obj, f.obj, rbind(f.con2, f.con.vrs), c(f.dir, "="), c(f.rhs2, 1))$objval
        
        # Intertemporal scores
        
        f.con12 <- cbind(f.con2[, 1], f.con1[, -1])
        f.con21 <- cbind(f.con1[, 1], f.con2[, -1])
        
        eff12[t - 1, i] <- lp(obj, f.obj, f.con12, f.dir, f.rhs2)$objval
        eff21[t - 1, i] <- lp(obj, f.obj, f.con21, f.dir, f.rhs1)$objval
        
        f.con1 <- f.con2
        f.rhs1 <- f.rhs2
        
      }
      
    }
    
  } else {
    
    obj <- "max"
    f.obj <- c(1, rep(0, ndr))
    
    f.dir <- c(rep(">=", no), rep("<=", ni))
    f.con.vrs <- cbind(0, matrix(1, nrow = 1, ncol = ndr))
    
    for (i in 1:nde) {
      
      ii <- dmu_eval[i]
      
      f.con.1 <- cbind(-datadealist[[1]]$output[, ii], 
                       matrix(datadealist[[1]]$output[, dmu_ref], nrow = no))
      f.con.2 <- cbind(matrix(0, nrow = ni, ncol = 1),
                       matrix(datadealist[[1]]$input[, dmu_ref], nrow = ni))
      f.con1 <- rbind(f.con.1, f.con.2)
      
      f.rhs1 <- c(rep(0, no), datadealist[[1]]$input[, ii])
      
      eff[1, i] <- lp(obj, f.obj, f.con1, f.dir, f.rhs1)$objval
      effvrs[1, i] <- lp(obj, f.obj, rbind(f.con1, f.con.vrs), c(f.dir, "="), c(f.rhs1, 1))$objval
      
      for (t in 2:nt) {
        
        f.con.1 <- cbind(-datadealist[[t]]$output[, ii],
                         matrix(datadealist[[t]]$output[, dmu_ref], nrow = no))
        f.con.2 <- cbind(matrix(0, nrow = ni, ncol = 1),
                         matrix(datadealist[[t]]$input[, dmu_ref], nrow = ni))
        f.con2 <- rbind(f.con.1, f.con.2)
        
        f.rhs2 <- c(rep(0, no), datadealist[[t]]$input[, ii])
        
        eff[t, i] <- lp(obj, f.obj, f.con2, f.dir, f.rhs2)$objval
        effvrs[t, i] <- lp(obj, f.obj, rbind(f.con2, f.con.vrs), c(f.dir, "="), c(f.rhs2, 1))$objval
        
        # Intertemporal scores
        
        f.con12 <- cbind(f.con2[, 1], f.con1[, -1])
        f.con21 <- cbind(f.con1[, 1], f.con2[, -1])
        
        eff12[t - 1, i] <- lp(obj, f.obj, f.con12, f.dir, f.rhs2)$objval
        eff21[t - 1, i] <- lp(obj, f.obj, f.con21, f.dir, f.rhs1)$objval
        
        f.con1 <- f.con2
        f.rhs1 <- f.rhs2
        
      }
      
    }
    
    eff <- 1 / eff
    effvrs <- 1 / effvrs
    eff12 <- 1 / eff12
    eff21 <- 1 / eff21
    
  }
  
  mi <- sqrt((eff12 * eff[-1, ]) / (eff[-nt, ] * eff21))
  ec <- eff[-1, ] / eff[-nt, ]
  tc <- mi / ec
  pech <- effvrs[-1, ]/effvrs[-nt, ]
  sech <- ec / pech
  
  deaOutput <- list(mi = mi,
              ec = ec,
              tc = tc,
              pech = pech,
              sech = sech,
              eff = eff,
              eff12 = eff12,
              eff21 = eff21,
              effvrs = effvrs,
              datadealist = datadealist,
              dmu_eval = dmu_eval,
              dmu_ref = dmu_ref,
              orientation = orientation,
              modelname = "malmquist")
  return(structure(deaOutput, class = "dea"))
  
}