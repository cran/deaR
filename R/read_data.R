#' @title read_data
#'  
#' @description This function creates, from a data frame, a \code{deadata} structure, which is as list with fields
#' \code{input}, \code{output}, \code{dmunames}, \code{nc_inputs}, \code{nc_outputs}, \code{nd_inputs}, \code{nd_outputs}.
#'
#' @usage read_data(datadea = NULL,
#'           ni = NULL,
#'           no = NULL,
#'           dmus = 1,
#'           inputs = NULL,
#'           outputs = NULL,
#'           nc_inputs = NULL,
#'           nc_outputs = NULL,
#'           nd_inputs = NULL,
#'           nd_outputs = NULL,
#'           ud_inputs = NULL, 
#'           ud_outputs = NULL)
#'              
#' @param datadea Data frame with DEA data.
#' @param dmus Column (number or name) of DMUs (optional). By default, it is the first column.
#' If there is not any DMU column, then it must be \code{NULL}.
#' @param ni Number of inputs, if inputs are in columns 2:(\code{ni} + 1) (if DMUs are in the first column) or 1:\code{ni} (no DMUs column).
#' @param no Number of outputs, if outputs are in columns (\code{ni} + 2):(\code{ni} + \code{no} + 1) (if DMUs are in the first column) or (\code{ni} + 1):(\code{ni} + \code{no}) (no DMUs column).
#' If not specified, DMUs are in the first column.
#' @param inputs Columns (numbers or names) of inputs (optional). It prevails over \code{ni}.
#' Alternatively to \code{datadea}, it can be a matrix with the inputs (DMUs in columns). In this case, DMUs names
#' are taken from the columns names.
#' @param outputs Columns (numbers or names) of outputs (optional). It prevails over \code{no}.
#' Alternatively to \code{datadea}, it can be a matrix with the outputs (DMUs in columns).
#' @param nc_inputs A numeric vector containing the indices of non-controllable inputs.
#' @param nc_outputs A numeric vector containing the indices of non-controllable outputs.
#' @param nd_inputs A numeric vector containing the indices of non-discretionary inputs.
#' @param nd_outputs A numeric vector containing the indices of non-discretionary outputs.
#' @param ud_inputs A numeric vector containing the indices of undesirable (good) inputs.
#' @param ud_outputs A numeric vector containing the indices of undesirable (bad) outputs.
#'
#' @return An object of class \code{deadata}
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
#' @examples
#' data("Coll_Blasco_2006")
#' data_example <- read_data(datadea = Coll_Blasco_2006,
#'                           ni = 2, 
#'                           no = 2)
#' # This is the same as:
#' data_example <- read_data(Coll_Blasco_2006,
#'                           inputs = 2:3, 
#'                           outputs = 4:5)
#' # And the same as:
#' dmunames <- c("A", "B", "C", "D", "E", "F")
#' nd <- length(dmunames) # Number of DMUs
#' inputnames <- c("Employees", "Capital")
#' ni <- length(inputnames) # Number of Inputs
#' outputnames <- c("Vehicles", "Orders")
#' no <- length(outputnames) # Number of Outputs
#' inputs <- matrix(c(8, 8, 11, 15, 14, 12, 12, 13, 11, 18, 18, 20),
#'                  nrow = ni, ncol = nd, dimnames = list(inputnames, dmunames))
#' outputs <- matrix(c(14, 20, 25, 42, 8, 30, 25, 8, 40, 22, 24, 30),
#'                   nrow = no, ncol = nd, dimnames = list(outputnames, dmunames))
#' data_example <- read_data(inputs = inputs,
#'                           outputs = outputs)                
#' # If the first input is a non-controllable input:
#' data_example <- read_data(Coll_Blasco_2006,
#'                           inputs = 2:3,
#'                           outputs = 4:5, 
#'                           nc_inputs = 1)
#' # If the second output is a non-discretionary output:
#' data_example <- read_data(Coll_Blasco_2006,
#'                           inputs = 2:3, 
#'                           outputs = 4:5, 
#'                           nd_outputs = 2)
#' # If the second input is a non-discretionary input and the second output is an undesirable:
#' data_example <- read_data(Coll_Blasco_2006,
#'                           inputs = 2:3, 
#'                           outputs = 4:5, 
#'                           nd_inputs = 2, 
#'                           ud_outputs = 2)
#' 
#' @export

read_data <- function(datadea = NULL,
                      ni = NULL,
                      no = NULL,
                      dmus = 1,
                      inputs = NULL,
                      outputs = NULL,
                      nc_inputs = NULL,
                      nc_outputs = NULL,
                      nd_inputs = NULL,
                      nd_outputs = NULL,
                      ud_inputs = NULL,
                      ud_outputs = NULL) {
  
  if (is.matrix(inputs) && is.matrix(outputs)) {
    nd <- ncol(inputs)
    if (ncol(outputs) != nd) {
      stop("Inputs and outputs matrices must have the same number of columns (number of DMUs).")
    }
    input <- inputs  
    output <- outputs
    dmunames <- colnames(inputs)
    if (is.null(dmunames)) {
      dmunames <- paste("DMU", 1:nd, sep = "")
      colnames(input) <- dmunames
    }
    colnames(output) <- dmunames
    ni <- nrow(inputs)
    no <- nrow(outputs)
    inputnames <- rownames(inputs)
    if (is.null(inputnames)) {
      inputnames <- paste("Input", 1:ni, sep = "")
      rownames(input) <- inputnames
    }
    outputnames <- rownames(outputs)
    if (is.null(outputnames)) {
      outputnames <- paste("Output", 1:no, sep = "")
      rownames(output) <- outputnames
    }
  } else {
    
    # Checking datadea
    if(!is.data.frame(datadea)){
      stop("Invalid data datadea (should be a data frame)!")
    }
    datadea <- as.data.frame(datadea)
    
    # Checking inputs
    if (is.null(inputs)){
      if (!is.null(ni)) {
        if (is.null(dmus)) {
          inputs <- 1:ni
        } else if (dmus == 1) {
          inputs <- 2:(ni + 1)
        } else {
          stop("If you specify ni, then dmus must be NULL or 1.")
        }
      } else {
        stop("Inputs not specified.")
      }
    }
    
    # Checking outputs
    if (is.null(outputs)){
      if (!is.null(no)) {
        if (is.null(dmus)) {
          outputs <- (ni + 1):(ni + no)
        } else if (dmus == 1) {
          outputs <- (ni + 2):(ni + no + 1)
        } else {
          stop("If you specify no, then dmus must be NULL or 1.")
        }
      } else {
        stop("Outputs not specified.")
      }
    }
    
    # Checking dmu
    if (is.null(dmus)) {
      dmunames <- paste("DMU", 1:nrow(datadea), sep = "")
    } else {
      if (length(dmus) > 1) {
        stop("Invalid DMU names specification. Provide either a single
             column number or name.")
      } else {
        if (!class(dmus) %in% c("integer", "numeric", "character")) {
          stop("Invalid DMU names specification. Please give either the column
               number or name.")
        } else {
          if (is.character(dmus) && !(dmus %in% colnames(datadea))) {
            stop(" Invalid DMU names. Please either give the DMU column number or name.")
          }
          if ((is.numeric(dmus) || is.integer(dmus)) && (dmus > ncol(datadea) || dmus < 1)) {
            stop("Invalid DMU names specification. Give NULL or a valid column number!")
          }
        }
    }
      
      dmunames <- datadea[, dmus]
  }
    
    if (is.character(inputs)) {
      inputnames <- inputs
    } else {
      inputnames <- colnames(datadea)[inputs]
    }
    if (is.character(outputs)) {
      outputnames <- outputs
    } else {
      outputnames <- colnames(datadea)[outputs]
    }
  
    ni <- length(inputnames)
    no <- length(outputnames)
    
    input <- datadea[, inputs]
    output <- datadea[, outputs]
    input <- t(input)  
    output <- t(output)
    rownames(input) <- inputnames
    rownames(output) <- outputnames
    colnames(input) <- dmunames
    colnames(output) <- dmunames
  
  }
  
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
  
  # Checking orders of magnitude in data
  maxio <- max(max(input), max(output), na.rm = TRUE)
  minio <- min(min(input), min(output), na.rm = TRUE)
  if (minio > 0) {
    if (maxio / minio > 1e6) {
      warning("There are data with very different orders of magnitude. Try to re-define the units of measure or some linear problems may be ill-posed.")
    }
  } else {
    warning("There are negative or zero data.")
  }

  res <- list(
    input = input,
    output = output,
    dmunames = dmunames,
    nc_inputs = nc_inputs,
    nc_outputs = nc_outputs,
    nd_inputs = nd_inputs,
    nd_outputs = nd_outputs,
    ud_inputs = ud_inputs,
    ud_outputs = ud_outputs
    )
  return(structure(res, class = "deadata"))
}