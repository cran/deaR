#' @title Extreme efficient DMUs
#'   
#' @description Find a set of extreme efficient DMUs from a \code{deadata} object.
#' 
#' @usage extreme_efficient(datadea,
#'              dmu_ref = NULL,
#'              rts = c("crs", "vrs", "nirs", "ndrs"),
#'              tol = 1e-6)
#' 
#' @param datadea A \code{deadata} object with \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set,
#'                i.e. the cluster of DMUs from which we want to find an extreme efficient DMUs subset.
#'                If \code{NULL} (default), all DMUs are considered.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#'            "vrs" (variable), "nirs" (non-increasing) or "ndrs" (non-decreasing).
#' @param tol Numeric, a tolerance margin for checking efficiency. It is 1e-6 by default.
#' 
#' @returns A numeric vector representing an extreme efficient subset of DMUs.
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
#' Charnes, A.; Cooper, W.W.; Thrall, R.M. (1991) "A structure for classifying
#' and characterizing efficiency and inefficiency in data envelopment analysis",
#' Journal of Productivity Analisys, 2, 197–237.
#' 
#' @examples
#' data("PFT1981")
#' datadea <- make_deadata(PFT1981,
#'                         ni = 5,
#'                         no = 3)
#' # We find an extreme efficient subset from a cluster formed by the first 20 DMUs
#' result <- extreme_efficient(datadea = datadea,
#'                             dmu_ref = 1:20)
#' 
#' @import lpSolve
#' 
#' @export

extreme_efficient <- function(datadea,
                              dmu_ref = NULL,
                              rts = c("crs", "vrs", "nirs", "ndrs"),
                              tol = 1e-6) {
  
  # Cheking whether datadea is of class "deadata" or not...  
  if (!is.deadata(datadea)) {
    stop("Data should be of class deadata. Run make_deadata function first!")
  }
  
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
  dmunames <- datadea$dmunames
  nd <- length(datadea$dmunames) # number of dmus
  
  input <- datadea$input
  output <- datadea$output
  ni <- nrow(input) # number of  inputs
  no <- nrow(output) # number of outputs
  
  if (is.null(dmu_ref)) {
    dmu_ref <- 1:nd
  } else if (!all(dmu_ref %in% (1:nd))) {
    stop("Invalid set of reference DMUs (dmu_ref).")
  }
  names(dmu_ref) <- dmunames[dmu_ref]
  
  # Find efficient DMUs in dmu_ref
  result_add <- model_additive(datadea = datadea,
                               dmu_eval = dmu_ref,
                               dmu_ref = dmu_ref,
                               rts = rts)
  # objval <- unlist(lapply(result_add$DMU, function(x) x$objval))
  slacksio <- slacks(result_add)
  slacks_input <- t(slacksio$slack_input) / datadea$input[, dmu_ref]
  slacks_output <- t(slacksio$slack_output) / datadea$output[, dmu_ref]
  slacks_matrix <- rbind(slacks_input, slacks_output)
  objval <- colSums(slacks_matrix)
  effDMUs <- dmu_ref[which(objval <= tol)]
  ne <- length(effDMUs)
  
  input_eff <- matrix(input[, effDMUs], nrow = ni)
  output_eff <- matrix(output[, effDMUs], nrow = no)
  
  extreme <- NULL
  f.obj <- rep(1, ne - 1)
  if (rts == "crs") {
    f.con.rs <- NULL
    f.dir.rs <- NULL
    f.rhs.rs <- NULL
  } else {
    f.con.rs <- f.obj
    f.rhs.rs <- 1
    if (rts == "vrs") {
      f.dir.rs <- "="
    } else if (rts == "nirs") {
      f.dir.rs <- "<="
    } else {
      f.dir.rs <- ">="
    }
  }
  f.dir <- c(rep("=", ni + no), f.dir.rs)
  # Si hay DMUs iguales o proporcionales hay que quedarse con una.
  # Lo hago al reves para ir quitando las ultimas y quedarme con la primera.
  for (i in ne:1) {
    ii <- effDMUs[i]
    input_i <- matrix(input[, effDMUs[-i]], nrow = ni)
    output_i <- matrix(output[, effDMUs[-i]], nrow = no)
    f.con <- rbind(input_i, output_i, f.con.rs)
    f.rhs <- c(input[, ii], output[, ii], f.rhs.rs)
    res <- lp("max", f.obj, f.con, f.dir, f.rhs)
    if (res$status != 0) {
      extreme <- c(ii, extreme)
    } else {
      effDMUs <- effDMUs[-i]
      f.obj <- f.obj[-1]
    }
  }
  
  return(extreme)
  
}