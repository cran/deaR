#' @title Range directional model.
#' 
#' @description Range directional model from Portela et al. (2004).
#' 
#' @note Undesirable inputs/outputs are treated as negative inputs/outputs in this model.
#' 
#' @usage model_rdm(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             orientation = c("no", "io", "oo"),
#'             irdm = FALSE,
#'             maxslack = TRUE,
#'             weight_slack_i = 1,
#'             weight_slack_o = 1,
#'             compute_target = TRUE,
#'             returnlp = FALSE,
#'             ...)
#' 
#' @param datadea A \code{deadata} object, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation A string, equal to "no" (non-oriented), "io" (input oriented),
#' or "oo" (output oriented).
#' @param irdm Logical. If it is \code{TRUE}, it applies the IRDM (inverse range
#' directional model).
#' @param maxslack Logical. If it is \code{TRUE}, it computes the max slack solution.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval}) with the weights
#' of the input slacks for the max slack solution.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x
#' \code{ne} (where \code{ne} is the length of \code{dmu_eval}) with the weights
#' of the output slacks for the max slack solution.
#' @param compute_target Logical. If it is \code{TRUE}, it computes targets of the
#' max slack solution. 
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
#' Portela, M.; Thanassoulis, E.; Simpson, G. (2004). "Negative data in DEA: a
#' directional distance approach applied to bank branches", Journal of the Operational
#' Research Society, 55 1111-1121.  
#' 
#' @export

model_rdm <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           orientation = c("no", "io", "oo"),
           irdm = FALSE,
           maxslack = TRUE,
           weight_slack_i = 1,
           weight_slack_o = 1,
           compute_target = TRUE,
           returnlp = FALSE,
           ...) {
    
    # Cheking whether datadea is of class "deadata" or not...  
    if (!is.deadata(datadea)) {
      stop("Data should be of class deadata. Run make_deadata function first!")
    }
      
    # Checking orientation
    orientation <- tolower(orientation)
    orientation <- match.arg(orientation)
    
    nd <- length(datadea$dmunames) # number of dmus
    
    if (is.null(dmu_eval)) {
      dmu_eval <- 1:nd
    } else if (!all(dmu_eval %in% (1:nd))) {
      stop("Invalid set of DMUs to be evaluated (dmu_eval).")
    }
    nde <- length(dmu_eval)
    
    if (is.null(dmu_ref)) {
      dmu_ref <- 1:nd
    } else if (!all(dmu_ref %in% (1:nd))) {
      stop("Invalid set of reference DMUs (dmu_ref).")
    }
    ndr <- length(dmu_ref)
    
    input <- datadea$input
    output <- datadea$output
    ni <- nrow(input) # number of  inputs
    no <- nrow(output) # number of outputs
    
    # Undesirable i/o
    ud_inputs <- datadea$ud_inputs
    ud_outputs <- datadea$ud_outputs
    if (!is.null(ud_inputs) || !is.null(ud_outputs)) {
      input[ud_inputs, ] <- -abs(input[ud_inputs, ])
      output[ud_outputs, ] <- -abs(output[ud_outputs, ])
    }
    
    # Directional vectors
    dir_input <- matrix(0, nrow = ni, ncol = nde)
    dir_output <- matrix(0, nrow = no, ncol = nde)
    
    if (orientation != "oo") {
      dir_input <- input[, dmu_eval] - apply(matrix(input[, dmu_ref], nrow = ni), MARGIN = 1, FUN = min)
      dir_input <- matrix(dir_input, nrow = ni)
    }
    if (orientation != "io") {
      dir_output <- apply(matrix(output[, dmu_ref], nrow = no), MARGIN = 1, FUN = max) - output[, dmu_eval]
      dir_output <- matrix(dir_output, nrow = no)
    }
    
    if (irdm) {
      dir_input[which(dir_input != 0)] <- 1 / dir_input[which(dir_input != 0)]
      dir_output[which(dir_output != 0)] <- 1 / dir_output[which(dir_output != 0)]
    }
    
    datadea_rdm <- datadea
    datadea_rdm$input <- input
    datadea_rdm$output <- output
    datadea_rdm$ud_inputs <- NULL
    datadea_rdm$ud_outputs <- NULL
    rdm <- model_basic(datadea = datadea_rdm, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                       orientation = "dir", dir_input = dir_input, dir_output = dir_output,
                       rts = "vrs", maxslack = maxslack, weight_slack_i = weight_slack_i, weight_slack_o = weight_slack_o,
                       compute_target = compute_target, returnlp = returnlp)
    
    if (!is.null(ud_inputs) || !is.null(ud_outputs)) {
      datadea_rdm$ud_inputs <- ud_inputs
      datadea_rdm$ud_outputs <- ud_outputs
    }
    
    deaOutput <- list(modelname = "rdm",
                     orientation = orientation,
                     rts = "vrs",
                     DMU = rdm$DMU,
                     data = datadea_rdm,
                     dmu_eval = rdm$dmu_eval,
                     dmu_ref = rdm$dmu_ref,
                     dir_input = rdm$orientation_param$dir_input,
                     dir_output = rdm$orientation_param$dir_output,
                     irdm = irdm,
                     maxslack = maxslack,
                     weight_slack_i = rdm$weight_slack_i,
                     weight_slack_o = rdm$weight_slack_o)
    
    return(structure(deaOutput, class = "dea"))
  }