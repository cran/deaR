#' @title Fuzzy DEA model
#'
#' @description Solve the Fuzzy input-oriented and output-oriented DEA model proposed by Guo and Tanaka (2001) under constant returns-to-scale.
#' In deaR is implemented the LP poblem given by the model (16) in Guo and Tanaka (2001, p.155).
#' The fuzzy efficiencies are calculated according to equations in (17) (Guo and Tanaka, 2001, p.155).
#' The (crisp) relative efficiencies and multipliers for the case \code{h} = 1 are obtained from the CCR model (\code{model_multiplier}).
#' @note The optimal solution of model (16) is not unique. 
#' @usage modelfuzzy_guotanaka(datadea,
#'                      dmu_eval = NULL,
#'                      dmu_ref = NULL,
#'                      orientation = c("io", "oo"),
#'                      h = 1)
#' 
#' @param datadea The data, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation A string, equal to "io" (input oriented) or "oo" (output oriented).
#' @param h A numeric vector with the h-levels (in [0,1]).
#'   
#' @return An object of class \code{deadata_fuzzy}.
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
#' Emrouznejad, A.; Tavana, M.; Hatami-Marbini, A. (2014). “The State of the Art in Fuzzy Data Envelopment Analysis”, in A. Emrouznejad and M. Tavana (eds.), Performance Measurement with Fuzzy Data Envelopment Analysis. Studies in Fuzziness and Soft Computing 309. Springer, Berlin. \doi{10.1007/978-3-642-41372-8_1}
#' 
#' Guo, P.; Tanaka, H. (2001). "Fuzzy DEA: A Perceptual Evaluation Method", Fuzzy Sets and Systems, 119, 149–160. \doi{10.1016/S0165-0114(99)00106-2}
#' 
#' Hatami-Marbini, A.; Emrouznejad, A.; Tavana, M. (2011). "A Taxonomy and Review of the Fuzzy Data Envelopment Analysis Literature: Two Decades in the Making", European Journal of Operational Research, 214, 457–472. \doi{10.1016/j.ejor.2011.02.001}
#' 
#' @examples
#' # Example 1.
#' # Replication results in Guo and Tanaka (2001, p. 159). 
#' # In deaR is implemented the LP poblem given by the model 16 in Guo and Tanaka (2001, p. 155).
#' # The fuzzy efficiencies are calculated according to equations in (17) (Guo and Tanaka, 2001,p.155).
#' data("Guo_Tanaka_2001")
#' data_example <- read_data_fuzzy(Guo_Tanaka_2001, 
#'                                 inputs.mL = 2:3, 
#'                                 inputs.dL = 4:5, 
#'                                 outputs.mL = 6:7,
#'                                 outputs.dL = 8:9)
#' result <- modelfuzzy_guotanaka(data_example, 
#'                                h = c(0, 0.5, 0.75, 1), 
#'                                orientation = "io")
#' efficiencies(result)
#'  
#' # Example 2. 
#' data("Guo_Tanaka_2001")
#' data_example <- read_data_fuzzy(Guo_Tanaka_2001, 
#'                                 inputs.mL = 2:3, 
#'                                 inputs.dL = 4:5, 
#'                                 outputs.mL = 6:7, 
#'                                 outputs.dL = 8:9)
#' result2 <- modelfuzzy_guotanaka(data_example, 
#'                                 h = seq(0, 1, by = 0.1), 
#'                                 orientation = "io")
#' efficiencies(result2)
#' 
#' @seealso \code{\link{model_basic}}, \code{\link{model_multiplier}}, \code{\link{modelfuzzy_kaoliu}}, \code{\link{modelfuzzy_possibilistic}}, \code{\link{cross_efficiency_fuzzy}}
#' 
#' @import lpSolve
#' 
#' @export
  
modelfuzzy_guotanaka <-
function(datadea,
         dmu_eval = NULL,
         dmu_ref = NULL,
         orientation = c("io", "oo"),
         h = 1) {
 
  # Cheking whether datadea is of class "deadata_fuzzy" or not...  
  if (!is.deadata_fuzzy(datadea)) {
    stop("Data should be of class deadata_fuzzy. Run read_data_fuzzy function first!")
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
    
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
    
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
  
  # Checking whether data is symmetrical triangular fuzzy
  if (any(datadea$input$mL != datadea$input$mR) || any(datadea$output$mL != datadea$output$mR) ||
      any(datadea$input$dL != datadea$input$dR) || any(datadea$output$dL != datadea$output$dR)){
    stop("Fuzzy data must be symmetric triangular.")
  }
  
  if (orientation == "io") {
    input.m <- datadea$input$mL # center
    output.m <- datadea$output$mL
    input.d <- datadea$input$dL # spread
    output.d <- datadea$output$dL
    orient <- 1
    objg0 <- "max"
    g0.dir <- c("=", "<=")
    f.dir <- c(">=", g0.dir, rep("<=", ndr + ndr))
  } else {
    input.m <- -datadea$output$mL # center
    output.m <- -datadea$input$mL
    input.d <- -datadea$output$dL # spread
    output.d <- -datadea$input$dL
    orient <- -1
    objg0 <- "min"
    g0.dir <- c("=", ">=")
    f.dir <- c("<=", g0.dir, rep("<=", ndr + ndr))
  }
  
  inputnames <- rownames(input.m)
  outputnames <- rownames(output.m)
  ni <- length(inputnames) # number of inputs
  no <- length(outputnames) # number of outputs
  
  # Checking h
  if (any(h > 1) || any(h < 0)){
    stop("Invalid h vector.")
  }
  h <- sort(unique(h))
  nh <- length(h) # number of h-levels
  hlevel <- vector(mode = "list", length = nh)
  names(hlevel) <- as.character(h)
  
  e_spread <- max(input.d / input.m)
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  if (h[nh] == 1) {
    datadea_crisp <- structure(list(input = input.m,
                                    output = output.m,
                                    dmunames = dmunames),
                              class = "deadata")
    deasol_crisp <- do.call(model_multiplier, list(datadea = datadea_crisp,
                                              dmu_eval = dmu_eval,
                                              dmu_ref = dmu_ref,
                                              orientation = orientation))
    for (j in 1:nde) {
      efficiency <- c(0, deasol_crisp$DMU[[j]]$efficiency, 0)
      names(efficiency) <- c("dL", "m", "dR")
      DMU[[j]] <- list(efficiency = efficiency,
                       multiplier_input = deasol_crisp$DMU[[j]]$multiplier_input,
                       multiplier_output = deasol_crisp$DMU[[j]]$multiplier_output)
    }
    hlevel[[nh]] <- list(input = list(Lower = input.m, Upper = input.m),
                         output = list(Lower = output.m, Upper = output.m),
                         DMU = DMU)
    nh <- nh - 1
  }
  
  if (nh > 0) {
    for (i in 1:nh) {
      
      # h-levels
      
      a <- h[i]
      input.L <- input.m - input.d * (1 - a)
      input.U <- input.m + input.d * (1 - a)
      output.L <- output.m - output.d * (1 - a)
      output.U <- output.m + output.d * (1 - a)
      inputref.L <- matrix(input.L[, dmu_ref], nrow = ni)
      inputref.U <- matrix(input.U[, dmu_ref], nrow = ni)
      outputref.L <- matrix(output.L[, dmu_ref], nrow = no)
      outputref.U <- matrix(output.U[, dmu_ref], nrow = no)
      
      one.L <- 1 - (1 - a) * e_spread
      one.U <- 1 + (1 - a) * e_spread
      
      g0.rhs <- orient * c(one.L, one.U)
      
      f.con.3 <- cbind(-t(inputref.L), t(outputref.L))
      f.con.4 <- cbind(-t(inputref.U), t(outputref.U))
      
      for (j in 1:nde) {
        
        jj <- dmu_eval[j]
        
        # Compute g0
        g0.obj <- input.d[, jj]
        g0.con <- rbind(input.L[, jj], input.U[, jj])
        g0 <- lp(objg0, g0.obj, g0.con, g0.dir, g0.rhs)$objval
        
        # Vector de coeficientes de la función objetivo
        f.obj <- c(rep(0, ni), output.L[, jj])
        
        # Matriz técnica
        f.con.1 <- c(g0.obj, rep(0, no))
        f.con.2 <- cbind(g0.con, matrix(0, nrow = 2, ncol = no))
        f.con <- rbind(f.con.1, f.con.2, f.con.3, f.con.4)
        
        # Vector de términos independientes
        f.rhs <- c(g0, orient * c(one.L, one.U), rep(0, ndr + ndr))
        
        res <- lp("max", f.obj, f.con, f.dir, f.rhs)$solution
        
        multiplier_input <- res[1 : ni]
        multiplier_output <- res[(ni + 1) : (ni + no)]
        names(multiplier_input) <- inputnames
        names(multiplier_output) <- outputnames
        
        eff.m <- (multiplier_output %*% output.m[, jj]) / (multiplier_input %*% input.m[, jj])
        eff.dL <- eff.m - (multiplier_output %*% output.L[, jj]) / (multiplier_input %*% input.U[, jj])
        eff.dR <- (multiplier_output %*% output.U[, jj]) / (multiplier_input %*% input.L[, jj]) - eff.m
        
        efficiency <- c(eff.dL, eff.m, eff.dR)
        names(efficiency) <- c("dL", "m", "dR")
        
        if (orientation == "io") {
          DMU[[j]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_input, multiplier_output = multiplier_output)
        } else {
          DMU[[j]] <- list(efficiency = efficiency,
                           multiplier_input = multiplier_output, multiplier_output = multiplier_input)
        }
      }
      
      if (orientation == "io") {
        hlevel[[i]] <- list(input = list(Lower = input.L, Upper = input.U),
                            output = list(Lower = output.L, Upper = output.U),
                            DMU = DMU)
      } else {
        hlevel[[i]] <- list(input = list(Lower = -output.L, Upper = -output.U),
                            output = list(Lower = -input.L, Upper = -input.U),
                            DMU = DMU)
      }
      
    }
  }
  
  deaOutput <- list(modelname = "fuzzy_guotanaka",
                    orientation = orientation,
                    h = h,
                    hlevel = hlevel,
                    data = datadea,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref)

  return(structure(deaOutput, class = "dea_fuzzy"))
  
}