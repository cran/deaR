#' @title Radial super-efficiency basic DEA model
#'   
#' @description Solve Andersen and Petersen radial Super-efficiency DEA model. 
#' 
#' @usage model_supereff(datadea,
#'                dmu_eval = NULL,
#'                dmu_ref = NULL,
#'                supereff_modelname = c("basic"),
#'                ...)
#' 
#' @param datadea An object of class \code{deadata}.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param supereff_modelname A string containing the name of the radial model to
#' apply super-efficiency.
#' @param ... \code{orientation}, \code{rts} and other model parameters.
#' 
#' @returns A list of class \code{dea} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{orientation}, \code{rts},
#'  \code{dmu_eval} and \code{dmu_ref}.
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
#' Andersen, P.; Petersen, N.C. (1993). "A procedure for ranking efficient units in
#' data envelopment analysis", Management Science, 39, 1261-1264. 
#' 
#' Tone, K. (2002). "A slacks-based measure of super-efficiency in data envelopment
#' analysis", European Journal of Operational Research, 143, 32-41.
#' 
#' @examples
#' # Example 1.
#' # Replication of results in Tone (2002, p.38)
#' data("Power_plants")
#' data_example <- make_deadata(Power_plants, 
#'                              ni = 4, 
#'                              no = 2)
#' result <- model_supereff(data_example, 
#'                          orientation = "io", 
#'                          rts = "crs") 
#' eff <- efficiencies(result)
#'  
#' # Example 2. 
#' # Results of Super-efficiency with vrs returns to scale show infeasibility solutions 
#' # for DMUs D4 and D6 (these DMUs are not shown in deaR results).
#' data("Power_plants")
#' data_example2 <- make_deadata(Power_plants, 
#'                               ni = 4, 
#'                               no = 2) 
#' result2 <- model_supereff(data_example2, 
#'                           orientation = "io", 
#'                           rts = "vrs") 
#' eff2 <- efficiencies(result2)
#'
#' @note (1) Radial super-efficiency model under variable (vrs, nirs, ndrs, grs)
#' returns to scale can be infeasible for certain DMUs. See example 2.
#' 
#' (2) DMUs with infeasible solution are not shown in the results.
#' 
#' @seealso \code{\link{model_basic}}, \code{\link{model_sbmsupereff}},
#' \code{\link{model_addsupereff}}
#' 
#' @export
  
model_supereff <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           supereff_modelname = c("basic"), #"deaps", "fdh", "multiplier"), This super-efficiency model also gives feasible solutions for these models.
           ...) {
    
  # Checking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
    
  optlist <- list(...)
  if ("orientation" %in% names(optlist)) {
    if (optlist$orientation == "dir") {
      stop("Directional models not supported.")
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
  
  # Checking modelname
  supereff_modelname <- tolower(supereff_modelname)
  supereff_modelname <- match.arg(supereff_modelname)
  model_modelname <- paste("model", supereff_modelname, sep = "_")
  
  # Zeros in input and output data. Case 2 (Tone 2001)
  nzimin <- apply(datadea$input, MARGIN = 1, function(x) min(x[x > 0])) / 100
  nzomin <- apply(datadea$output, MARGIN = 1, function(x) min(x[x > 0])) / 100
  for (ii in dmu_eval) {
    datadea$input[which(datadea$input[, ii] == 0), ii] <- nzimin[which(datadea$input[, ii] == 0)]
    datadea$output[which(datadea$output[, ii] == 0), ii] <- nzomin[which(datadea$output[, ii] == 0)]
  }
  
  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]
  
  for (i in 1:nde) {
    
    ii <- dmu_eval[i]
    
    deasol <- do.call(model_modelname,
                      list(datadea = datadea,
                           dmu_eval = ii,
                           dmu_ref = dmu_ref[dmu_ref != ii],
                           ...))
    
    DMU[[i]] <- deasol$DMU[[1]]
    
    if ((ii %in% dmu_ref) && (!is.null(DMU[[i]]$lambda))) {
      newlambda <- rep(0, ndr)
      # newlambda[dmu_ref == ii] <- 0
      newlambda[dmu_ref != ii] <- DMU[[i]]$lambda
      names(newlambda) <- dmunames[dmu_ref]
      DMU[[i]]$lambda <- newlambda
    }
    
  }
  
  deaOutput <- deasol
  
  deaOutput$modelname <- paste("supereff", deasol$modelname, sep = "_")
  deaOutput$DMU <- DMU
  deaOutput$dmu_eval <- dmu_eval
  deaOutput$dmu_ref <- dmu_ref
 
  return(deaOutput)
 
}
