#' @title Free disposal hull (FDH) model.
#'   
#' @description  FDH model allows the free disposability to construct the production
#' possibility set. The central feature of the FDH model is the lack of convexity
#' for its production possibility set (Thrall, 1999).
#' 
#' @usage model_fdh(datadea,
#'           fdh_modelname = c("basic"),
#'           ...)
#' 
#' @param datadea A \code{deadata} object, including DMUs, inputs and outputs.
#' @param fdh_modelname A string containing the name of the model to apply FDH.
#' For now, only "basic" is available.
#' @param ... \code{dmu_eval}, \code{dmu_ref}, \code{orientation} and other model parameters.
#' Parameters like \code{rts}, \code{max_slack} and \code{returnlp} are ignored.
#' 
#' @returns A list of class \code{dea} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model.
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
#' Cherchye, L.; Kuosmanen, T.; Post, T. (2000). "What Is the Economic Meaning of
#' FDH? A Reply to Thrall". Journal of Productivity Analysis, 13(3), 263–267. 
#' 
#' Deprins, D.; Simar, L. and Tulkens, H. (1984). Measuring Labor-Efficiency in
#' Post Offices. In M. Marchand, P. Pestieau and H. Tulkens (eds.), The Performance
#' of Public Entreprises: Concepts and Measurement. Amsterdam: North-Holland.
#' 
#' Sanei, M.; Mamizadeh Chatghayeb, S. (2013). “Free Disposal Hull Models in
#' Supply Chain Management”, International Journal of Mathematical Modelling and
#' Computations, 3(3), 125-129.
#' 
#' Thrall, R. M. (1999).  "What Is the Economic Meaning of FDH?", Journal of
#' Productivity Analysis, 11(3), 243–50. 
#' 
#' @examples 
#' # Example 1. FDH input-oriented.
#' # Replication of results in Sanei and Mamizadeh Chatghayeb (2013)
#' data("Supply_Chain")
#' data_fdh1 <- make_deadata(Supply_Chain, 
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- model_fdh(data_fdh1) # by default orientation = "io"
#' efficiencies(result)
#' 
#' # Example 2. FDH output-oriented.
#' # Replication of results in Sanei and Mamizadeh Chatghayeb (2013)
#' data("Supply_Chain")
#' data_fdh2 <- make_deadata(Supply_Chain, 
#'                           inputs = 5:6, 
#'                           outputs = 7:8)
#' result2 <- model_fdh(data_fdh2, 
#'                     orientation = "oo")
#' efficiencies(result2)
#'  
#' @import lpSolve
#' 
#' @export
  
model_fdh <-
  function(datadea,
           fdh_modelname = c("basic"), # It can be generalized to other models in the future
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
    
  # Checking modelname
  fdh_modelname <- tolower(fdh_modelname)
  fdh_modelname <- match.arg(fdh_modelname)
  model_modelname <- paste("model", fdh_modelname, sep = "_")
  optlist <- list(...)
  if ("rts" %in% names(optlist)) {
    #optlist$rts <- "vrs" # lpSolve and Rglpk do not work properly with binary variables
    optlist$rts <- "crs"
  }
  optlist$returnlp <- TRUE
  optlist$datadea <- datadea
  deasol <- do.call(model_modelname,
                    optlist)
  
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    if (fdh_modelname == "basic") {
      res_und <- undesirable_basic(datadea = datadea, vtrans_i = deasol$vtrans_i,
                                   vtrans_o = deasol$vtrans_o)
      datadea <- res_und$u_datadea
    }
  }
  
  input <- datadea$input
  output <- datadea$output
  dmunames <- datadea$dmunames
  inputnames <- rownames(input)
  outputnames <- rownames(output)
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  dmu_ref <- deasol$dmu_ref
  ndr <- length(dmu_ref)
  inputref <- matrix(input[, dmu_ref], nrow = ni) 
  outputref <- matrix(output[, dmu_ref], nrow = no)
  
  ncd_inputs <- c(datadea$nc_inputs, datadea$nd_inputs)
  ncd_outputs <- c(datadea$nc_outputs, datadea$nd_outputs)
  
  dmu_eval <- deasol$dmu_eval
  nde <- length(dmu_eval)
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  binary.vec = (2 : (ndr + 1)) # lambdas. This is valid for model_basic
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    obj <- deasol$DMU[[i]]$direction
    f.obj <- deasol$DMU[[i]]$objective.in
    f.con <- deasol$DMU[[i]]$const.mat
    f.dir <- deasol$DMU[[i]]$const.dir
    f.rhs <- deasol$DMU[[i]]$const.rhs
    
    # lpSolve and Rglpk do not work properly with binary variables
    # res <- lp(obj, f.obj, f.con, f.dir, f.rhs, binary.vec = binary.vec)$solution
    # kkmax <- (obj == "max")
    # res <- Rglpk_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, rhs = f.rhs,
    #                       types = c("C", rep("B", ndr)), max = kkmax)$solution
    
    ### Alternative ###
    
    scoreskk <- rep(NA, ndr)
    f.objkk <- f.obj[1]
    f.conkk <- matrix(f.con[, 1], ncol = 1)
    
    for (j in 1:ndr) {
      f.rhskk <- f.rhs - f.con[, j + 1]
      reskk <- lp(obj, f.objkk, f.conkk, f.dir, f.rhskk)
      if (reskk$status == 0) {
        scoreskk[j] <- reskk$solution[1]
      }
    }
    res <- rep(0, ndr + 1)
    if (obj == "min") {
      res[1] <- min(scoreskk, na.rm = TRUE)
    } else {
      res[1] <- max(scoreskk, na.rm = TRUE)
    }
    jkk <- which(scoreskk == res[1])[1]
    res[jkk + 1] <- 1
    
    ### End Alternative ###
    
    if (fdh_modelname == "basic") { # Prepared for generalization to other models
      
      lambda <- res[binary.vec]
      names(lambda) <- dmunames[dmu_ref]
      
      target_input <- as.vector(inputref %*% lambda)
      names(target_input) <- inputnames
      target_output <- as.vector(outputref %*% lambda)
      names(target_output) <- outputnames
      
      if (deasol$orientation == "io") {
        
        efficiency <- res[1]
        
        slack_input <- efficiency * input[, ii] - target_input
        slack_input[ncd_inputs] <- input[, ii] - target_input
        names(slack_input) <- inputnames
        slack_output <- target_output - output[, ii]
        names(slack_output) <- outputnames
        
        target_input[datadea$ud_inputs] <- deasol$vtrans_i - target_input[datadea$ud_inputs]
        target_output[datadea$ud_outputs] <- deasol$vtrans_o - target_output[datadea$ud_outputs]
        
        DMU[[i]] <- list(efficiency = efficiency,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output)
        
      } else if (deasol$orientation == "oo") {
        
        efficiency <- res[1]
        
        slack_input <- input[, ii] - target_input
        names(slack_input) <- inputnames
        slack_output <- target_output - efficiency * output[, ii]
        slack_output[ncd_outputs] <- target_output - output[, ii]
        names(slack_output) <- outputnames
        
        target_input[datadea$ud_inputs] <- deasol$vtrans_i - target_input[datadea$ud_inputs]
        target_output[datadea$ud_outputs] <- deasol$vtrans_o - target_output[datadea$ud_outputs]
        
        DMU[[i]] <- list(efficiency = efficiency,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output)
        
        
      } else {
        
        beta <- res[1]
        
        slack_input <- input[, ii] - beta * deasol$orientation_param$dir_input[, ii] - target_input
        slack_input[ncd_inputs] <- input[, ii] - target_input
        names(slack_input) <- inputnames
        slack_output <- target_output - output[, ii] - beta * deasol$orientation_param$dir_output[, ii]
        slack_output[ncd_outputs] <- target_output - output[, ii]
        names(slack_output) <- outputnames
        
        target_input[datadea$ud_inputs] <- deasol$vtrans_i - target_input[datadea$ud_inputs]
        target_output[datadea$ud_outputs] <- deasol$vtrans_o - target_output[datadea$ud_outputs]
        
        DMU[[i]] <- list(beta = beta,
                         lambda = lambda,
                         slack_input = slack_input, slack_output = slack_output,
                         target_input = target_input, target_output = target_output)
        
      }
      
    }
    
  }
  
  deaOutput <- deasol
  deaOutput$modelname <- paste("fdh", deasol$modelname, sep = "_")
  deaOutput$DMU <- DMU
 
  return(deaOutput)
 
}
