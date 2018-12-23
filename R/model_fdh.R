#' @title Free disposal hull (FDH) model.
#'   
#' @description  FDH model allows the free disposability to construct the production possiblity set. The central feature of the FDH model is the lack of convexity for its production possibility set (Thrall, 1999).
#' 
#' @usage model_fdh(datadea,
#'             fdh_modelname = c("basic"),
#'             ...)
#' 
#' @param datadea The data, including DMUs, inputs and outputs.
#' @param fdh_modelname A string containing the name of the model to apply FDH.
#' @param ... \code{dmu_eval}, \code{dmu_ref}, \code{orientation} and other model parameters.
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
#' Cherchye, L.; Kuosmanen, T.; Post, T. (2000). "What Is the Economic Meaning of FDH? A Reply to Thrall". Journal of Productivity Analysis, 13(3), 263–267. 
#' 
#' Deprins, D.; Simar, L. and Tulkens, H. (1984). Measuring Labor-Efficiency in Post Offices. In M. Marchand, P. Pestieau and H. Tulkens (eds.), The Performance of Public Entreprises: Concepts and Measurement. Amsterdam: North-Holland.
#' 
#' Sanei, M.; Mamizadeh Chatghayeb, S. (2013). “Free Disposal Hull Models in Supply Chain Management”, International Journal of Mathematical Modelling and Computations, 3(3), 125-129.
#' 
#' Thrall, R. M. (1999).  "What Is the Economic Meaning of FDH?", Journal of Productivity Analysis, 11(3), 243–50. 
#' 
#' @examples 
#' # Example 1. FDH input-oriented.
#' # Replication of results in Sanei and Mamizadeh Chatghayeb (2013)
#' data("Supply_Chain")
#' data_fdh1 <- read_data(Supply_Chain, 
#'                        dmus = 1, 
#'                        inputs = 2:4, 
#'                        outputs = 5:6)
#' result <- model_fdh(data_fdh1) # by default orientation = "io"
#' efficiencies(result)
#' 
#' # Example 2. FDH output-oriented.
#' # Replication of results in Sanei and Mamizadeh Chatghayeb (2013)
#' data("Supply_Chain")
#' data_fdh2 <- read_data(Supply_Chain, 
#'                        dmus = 1, 
#'                        inputs = 5:6, 
#'                        outputs = 7:8)
#' result2 <- model_fdh(data_fdh2, 
#'                     orientation = "oo")
#' efficiencies(result2)
#'  
#' @import lpSolve
#' 
#' @export
  
model_fdh <-
  function(datadea,
           fdh_modelname = c("basic"),
           ...) {
    
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run read_data function first!")
  }
    
  # Checking modelname
  fdh_modelname <- tolower(fdh_modelname)
  fdh_modelname <- match.arg(fdh_modelname)
  model_modelname <- paste("model", fdh_modelname, sep = "_")
    
  deasol <- do.call(model_modelname,
                    list(datadea = datadea,
                         rts = "vrs",
                         returnlp = TRUE,
                         ...))
  
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    if (fdh_modelname == "basic") {
      res_und <- undesirable_basic(datadea = datadea, vtrans_i = deasol$vtrans_i, vtrans_o = deasol$vtrans_o)
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
  
  binary.vec = (2 : (ndr + 1))
  
  dmu_eval <- deasol$dmu_eval
  nde <- length(dmu_eval)
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    obj <- deasol$DMU[[i]]$direction
    f.obj <- deasol$DMU[[i]]$objective.in
    f.con <- deasol$DMU[[i]]$const.mat
    f.dir <- deasol$DMU[[i]]$const.dir
    f.rhs <- deasol$DMU[[i]]$const.rhs
    
    res <- lp(obj, f.obj, f.con, f.dir, f.rhs, binary.vec = binary.vec)$solution
    
    if (fdh_modelname == "basic") {
      
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
        
        slack_input <- input[, ii] - beta * deasol$orientation$dir_input[, ii] - target_input
        slack_input[ncd_inputs] <- input[, ii] - target_input
        names(slack_input) <- inputnames
        slack_output <- target_output - output[, ii] - beta * deasol$orientation$dir_output[, ii]
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
