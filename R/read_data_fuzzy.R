#' @title read_data_fuzzy
#'  
#' @description This function creates, from a data frame, a \code{deadata_fuzzy} structure, which is as list with fields
#'   \code{input}, \code{output} and \code{dmunames}. At the same time, \code{input} and \code{output} are lists with fields
#'   \code{mL}, \code{mR}, \code{dL} and \code{dR}.
#'   \figure{fuzzynumbers.jpg}{options: width="100\%" alt="Figure: fuzzynumbers.jpg"}

#' @usage read_data_fuzzy(datadea,
#'              dmus,
#'              inputs.mL,
#'              inputs.mR = NULL,
#'              inputs.dL = NULL,
#'              inputs.dR = NULL,
#'              outputs.mL,
#'              outputs.mR = NULL,
#'              outputs.dL = NULL,
#'              outputs.dR = NULL,
#'              nc_inputs = NULL,
#'              nc_outputs = NULL,
#'              nd_inputs = NULL,
#'              nd_outputs = NULL,
#'              ud_inputs = NULL,
#'              ud_outputs = NULL)
#'          
#' @param datadea Data frame, the data.
#' @param dmus Where are (column) the DMUs in \code{datadea}.
#' @param inputs.mL Where are (columns) the \code{mL} of fuzzy inputs in \code{datadea}. If an input is not fuzzy or triangular, we put the value or \code{m} respectively.
#' @param inputs.mR Where are (columns) the \code{mR} of fuzzy inputs in \code{datadea}. If an input is not fuzzy or triangular, we put \code{NA}.
#' @param inputs.dL Where are (columns) the \code{dL} of fuzzy inputs in \code{datadea}. If an input is symmetric, we put \code{d}. If an input is not fuzzy or rectangular, we put \code{NA}.
#' @param inputs.dR Where are (columns) the \code{dR} of fuzzy inputs in \code{datadea}. If an input is not fuzzy or symmetric, we put \code{NA}.
#' @param outputs.mL Where are (columns) the \code{mL} of fuzzy outputs in \code{datadea}. If an output is not fuzzy or triangular, we put the value or \code{m} respectively.
#' @param outputs.mR Where are (columns) the \code{mR} of fuzzy outputs in \code{datadea}. If an output is not fuzzy or triangular, we put \code{NA}.
#' @param outputs.dL Where are (columns) the \code{dL} of fuzzy outputs in \code{datadea}. If an output is symmetric, we put \code{d}. If an output is not fuzzy or rectangular, we put \code{NA}.
#' @param outputs.dR Where are (columns) the \code{dR} of fuzzy outputs in \code{datadea}. If an output is not fuzzy or symmetric, we put \code{NA}.
#' @param nc_inputs A numeric vector containing non-controllable inputs.
#' @param nc_outputs A numeric vector containing non-controllable outputs.
#' @param nd_inputs A numeric vector containing non-discretionary inputs.
#' @param nd_outputs A numeric vector containing non-discretionary outputs.
#' @param ud_inputs A numeric vector containing undesirable (good) inputs.
#' @param ud_outputs A numeric vector containing undesirable (bad) outputs.
#'
#' @return An object of class \code{deadata_fuzzy}.
#'  
#' @examples
#' 
#' # Example 1. If inputs and/or outputs are symmetric triangular fuzzy numbers
#' data("Leon2003")
#' data_example <- read_data_fuzzy(datadea = Leon2003, 
#'                                 dmus = 1,
#'                                 inputs.mL = 2,
#'                                 inputs.dL = 3,
#'                                 outputs.mL = 4,
#'                                 outputs.dL = 5)
#' # Example 2. If inputs and/or outputs are non-symmetric triangular fuzzy numbers
#' data("Kao_Liu_2003")
#' data_example <- read_data_fuzzy(Kao_Liu_2003, 
#'                                 dmus = 1, 
#'                                 inputs.mL = 2, 
#'                                 outputs.mL = 3:7, 
#'                                 outputs.dL = c(NA, NA, 8, NA, 10),
#'                                 outputs.dR = c(NA, NA, 9, NA, 11))
#'                                 
#' @export

read_data_fuzzy <- function(datadea,
                            dmus,
                            inputs.mL,
                            inputs.mR = NULL,
                            inputs.dL = NULL,
                            inputs.dR = NULL,
                            outputs.mL,
                            outputs.mR = NULL,
                            outputs.dL = NULL,
                            outputs.dR = NULL,
                            nc_inputs = NULL,
                            nc_outputs = NULL,
                            nd_inputs = NULL,
                            nd_outputs = NULL,
                            ud_inputs = NULL,
                            ud_outputs = NULL) {
  
  # Checking data...
  if(!is.data.frame(datadea)){
    stop("Invalid input data (should be a data frame)!")
  }
  
  # Checking dmu
  if (is.null(dmus)) {
    dmunames <- paste("DMU", seq_along(datadea), sep = "")
  } else{
    if (length(dmus) > 1) {
      stop("Invalid dmu names specification. Provide either a single
           column number or name.")
    } else {
      if (!class(dmus) %in% c("integer", "numeric", "character")) {
        stop("Invalid dmu names specification. Please give either the column
             number or name.")
      } else{
        if (is.character(dmus) & !(dmus %in% colnames(datadea))) {
          stop(" Invalid dmu names. Please either give the dmu column number or name.")
        }
        if((is.numeric(dmus) | is.integer(dmus)) & (dmus > ncol(datadea) | dmus < 1)){
          stop("Invalid dmu names specification. Give a column number (dmus > ncols(datadea))!")
        }
      }
    }
    
    dmunames <- datadea[, dmus]
  }
  
  nd <- length(dmunames) # Number of DMUs
  ni <- length(inputs.mL) # Number of inputs
  no <- length(outputs.mL) # Number of outputs
  
  if(!is.character(inputs.mL)) {
    inputnames <- colnames(datadea)[inputs.mL]
  } else {
    inputnames <- inputs.mL
  }
  if(!is.character(outputs.mL)) {
    outputnames <- colnames(datadea)[outputs.mL]
  } else {
    outputnames <- outputs.mL
  }
  input.mL <- t(datadea[, inputs.mL])
  output.mL <- t(datadea[, outputs.mL])
  colnames(input.mL) <- dmunames
  colnames(output.mL) <- dmunames
  rownames(input.mL) <- inputnames
  rownames(output.mL) <- outputnames

  if (is.null(inputs.mR)) {
    input.mR <- input.mL
  } else {
    inputs.mR[is.na(inputs.mR)] <- inputs.mL[is.na(inputs.mR)] # NA in inputs.mR replaced by values of inputs.mL
    input.mR <- t(datadea[, inputs.mR])
  }
  colnames(input.mR) <- dmunames
  rownames(input.mR) <- inputnames

  if (is.null(outputs.mR)) {
    output.mR <- output.mL
  } else {
    outputs.mR[is.na(outputs.mR)] <- outputs.mL[is.na(outputs.mR)] # NA in outputs.mR replaced by values of outputs.mL
    output.mR <- t(datadea[, outputs.mR])
  }
  colnames(output.mR) <- dmunames
  rownames(output.mR) <- outputnames
  
  input.dL <- matrix(0, nrow = ni, ncol = nd)
  if (!is.null(inputs.dL)) {
    input.dL[!is.na(inputs.dL), ] <- t(datadea[, inputs.dL[!is.na(inputs.dL)]])
  }
  colnames(input.dL) <- dmunames
  rownames(input.dL) <- inputnames
  
  input.dR <- input.dL
  if (!is.null(inputs.dR)) {
    input.dR[!is.na(inputs.dR), ] <- t(datadea[, inputs.dR[!is.na(inputs.dR)]])
  }
  rownames(input.dR) <- inputnames
  
  output.dL <- matrix(0, nrow = no, ncol = nd)
  if (!is.null(outputs.dL)) {
    output.dL[!is.na(outputs.dL), ] <- t(datadea[, outputs.dL[!is.na(outputs.dL)]])
  }
  colnames(output.dL) <- dmunames
  rownames(output.dL) <- outputnames
  
  output.dR <- output.dL
  if (!is.null(outputs.dR)) {
    output.dR[!is.na(outputs.dR), ] <- t(datadea[, outputs.dR[!is.na(outputs.dR)]])
  }
  rownames(output.dR) <- outputnames
  
  # Checking non-controllable inputs/outputs
  if ((!is.null(nc_inputs)) && (!all(nc_inputs %in% 1:ni))) {
    stop("Invalid set of non-controllable inputs.")
  }
  if ((!is.null(nc_outputs)) && (!all(nc_outputs %in% 1:no))) {
    stop("Invalid set of non-controllable outputs.")
  }
  if (!is.null(nc_inputs)) {
    names(nc_inputs) <- inputnames[nc_inputs]
  }
  if (!is.null(nc_outputs)) {
    names(nc_outputs) <- outputnames[nc_outputs]
  }
  
  # Checking non-discretionary inputs/outputs
  if ((!is.null(nd_inputs)) && ((!all(nd_inputs %in% 1:ni)) || (any(nd_inputs %in% nc_inputs)))) {
    stop("Invalid set of non-discretionary inputs.")
  }
  if ((!is.null(nd_outputs)) && ((!all(nd_outputs %in% 1:no)) || (any(nd_outputs %in% nc_outputs)))) {
    stop("Invalid set of non-discretionary outputs.")
  }
  if (!is.null(nd_inputs)) {
    names(nd_inputs) <- inputnames[nd_inputs]
  }
  if (!is.null(nd_outputs)) {
    names(nd_outputs) <- outputnames[nd_outputs]
  }
  
  # Checking undesirable inputs/outputs
  if ((!is.null(ud_inputs)) && (!all(ud_inputs %in% 1:ni))) {
    stop("Invalid set of undesirable inputs.")
  }
  if ((!is.null(ud_outputs)) && (!all(ud_outputs %in% 1:no))) {
    stop("Invalid set of undesirable outputs.")
  }
  if (!is.null(ud_inputs)) {
    names(ud_inputs) <- inputnames[ud_inputs]
  }
  if (!is.null(ud_outputs)) {
    names(ud_outputs) <- outputnames[ud_outputs]
  }
  
  res <- list(
    input = list(
      mL = input.mL,
      mR = input.mR,
      dL = input.dL,
      dR = input.dR
    ),
    output = list(
      mL = output.mL,
      mR = output.mR,
      dL = output.dL,
      dR = output.dR
    ),
    dmunames = dmunames,
    nc_inputs = nc_inputs,
    nc_outputs = nc_outputs,
    nd_inputs = nd_inputs,
    nd_outputs = nd_outputs,
    ud_inputs = ud_inputs,
    ud_outputs = ud_outputs
  )
  
  return(structure(res, class = "deadata_fuzzy"))
}