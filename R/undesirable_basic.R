#' @title Undesirable inputs and outputs for basic DEA model.
#'  
#' @description This function transforms, according to Seiford and Zhu (2002),
#' a deadata or deadata_fuzzy class with undesirable inputs/outputs into a
#' deadata or deadata_fuzzy class with no undesirable inputs/outputs.
#' Onwards, it is recommended to use a DEA model with variable returns to scale (vrs).
#'
#' @usage undesirable_basic(datadea,
#'                   vtrans_i = NULL,
#'                   vtrans_o = NULL)
#' 
#' @param datadea The data, including DMUs, inputs and outputs.
#' @param vtrans_i Numeric vector of translation for undesirable inputs. If \code{vtrans_i[i]} is
#' \code{NA}, then it applies the "max + 1" translation to the i-th undesirable input.
#' If \code{vtrans_i} is a constant, then it applies the same translation to all
#' undesirable inputs. If \code{vtrans_i} is \code{NULL}, then it applies the
#' "max + 1" translation to all undesirable inputs.
#' @param vtrans_o Numeric vector of translation for undesirable outputs, analogous to
#'  \code{vtrans_i}, but applied to outputs.
#'
#' @return An object of class \code{deadata} or \code{deadata_fuzzy}.
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
#' Seiford, L.M.; Zhu, J. (2002). “Modeling undesirable factors in efficiency evaluation”,
#' European Journal of Operational Research 142, 16-20.
#' 
#' Hua Z.; Bian Y. (2007). DEA with Undesirable Factors. In: Zhu J., Cook W.D. (eds)
#' Modeling Data Irregularities and Structural Complexities in Data Envelopment Analysis.
#' Springer, Boston, MA. 
#' 
#' @examples
#' data("Hua_Bian_2007")
#' # The third output is an undesirable output.
#' data_example <- read_data(Hua_Bian_2007,
#'                           ni = 2, 
#'                           no = 3, 
#'                           ud_outputs = 3) 
#' # rts must be "vrs" for undesirable inputs/outputs:
#' # Translation parameter is set to (max + 1)
#' result <- model_basic(data_example,
#'                       orientation = "oo", 
#'                       rts = "vrs") 
#'
#' @export

undesirable_basic <- function(datadea,
                           vtrans_i = NULL,
                           vtrans_o = NULL) {
  
  if (is.deadata(datadea)) {
  
    # Undesirable Inputs
  
    input <- datadea$input
    ni <- nrow(input)
  
    ud_inputs <- datadea$ud_inputs
    nui <- length(ud_inputs)
  
    if (!all(ud_inputs %in% 1:ni)) {
      stop("Invalid set of undesirable inputs.")
    }
  
    if (is.null(vtrans_i)) {
      vtrans_i <- rep(NA, nui)
    } else if ((length(vtrans_i) == 1) && (nui > 1)) {
      vtrans_i <- rep(vtrans_i, nui)
    } else if (length(vtrans_i) != nui) {
      stop("Translation vector vtrans_i must be NULL, a constant or a vector of the
           same length as ud_inputs.")
    }
  
    if (nui > 0) {
      for (i in 1:nui) {
        if (is.na(vtrans_i[i])) {
          vtrans_i[i] <- max(input[ud_inputs[i], ]) + 1
        }
      }
    }
  
    i <- 1
    for (j in ud_inputs) {
      input[j, ] <- vtrans_i[i] - input[j, ]
      if (length(input[j, input[j, ] < 0] > 0)) {
        stop("There is at least one negative good input. Change translation vector vtrans_i.")
      }
      i <- i + 1
    }
  
    # Undesirable outputs
  
    output <- datadea$output
    no <- nrow(output)
  
    ud_outputs <- datadea$ud_outputs
    nuo <- length(ud_outputs)
  
    if (!all(ud_outputs %in% 1:no)) {
      stop("Invalid set of undesirable outputs.")
    }
  
    if (is.null(vtrans_o)) {
      vtrans_o <- rep(NA, nuo)
    } else if ((length(vtrans_o) == 1) && (nuo > 1)) {
      vtrans_o <- rep(vtrans_o, nuo)
    } else if (length(vtrans_o) != nuo) {
      stop("Translation vector vtrans_o must be NULL, a constant or a vector of the
           same length as ud_outputs.")
    }
  
    if (nuo > 0) {
      for (i in 1:nuo) {
        if (is.na(vtrans_o[i])) {
          vtrans_o[i] <- max(output[ud_outputs[i], ]) + 1
        }
      }
    }
  
    i <- 1
    for (j in ud_outputs) {
      output[j, ] <- vtrans_o[i] - output[j, ]
      if (length(output[j, output[j, ] < 0] > 0)) {
        stop("There is at least one negative bad output. Change translation vector vtrans_o.")
      }
      i <- i + 1
    }
    
    u_datadea <- datadea
    u_datadea$input <- input
    u_datadea$output <- output
    
    return(list(
      u_datadea = structure(u_datadea, class = "deadata"),
      vtrans_i = vtrans_i,
      vtrans_o = vtrans_o))
  
  } else if (is.deadata_fuzzy(datadea)) {
    
    # Undesirable Inputs
    
    input <- datadea$input
    ni <- nrow(input$mL)
    
    ud_inputs <- datadea$ud_inputs
    nui <- length(ud_inputs)
    
    if (!all(ud_inputs %in% 1:ni)) {
      stop("Invalid set of undesirable inputs.")
    }
    
    if (is.null(vtrans_i)) {
      vtrans_i <- rep(NA, nui)
    } else if ((length(vtrans_i) == 1) && (nui > 1)) {
      vtrans_i <- rep(vtrans_i, nui)
    } else if (length(vtrans_i) != nui) {
      stop("Translation vector vtrans_i must be NULL, a constant or a vector of the
           same length as ud_inputs.")
    }
    
    if (nui > 0) {
      for (i in 1:nui) {
        if (is.na(vtrans_i[i])) {
          vtrans_i[i] <- max(input$mR[ud_inputs[i], ] + input$dR[ud_inputs[i], ]) + 1
        }
      }
    }
    
    i <- 1
    for (j in ud_inputs) {
      
      aux <- input$dL[j, ]
      input$dL[j, ] <- input$dR[j, ]
      input$dR[j, ] <- aux
      
      aux <- input$mR[j, ]
      input$mR[j, ] <- vtrans_i[i] - input$mL[j, ]
      input$mL[j, ] <- vtrans_i[i] - aux
      
      if (length(input$mL[j, (input$mL[j, ] - input$dL[j, ]) < 0] > 0)) {
        stop("There is at least one negative good input. Change translation vector vtrans_i.")
        print(j)
      }
      
      i <- i + 1
    }
    
    # Undesirable outputs
    
    output <- datadea$output
    no <- nrow(output$mL)
    
    ud_outputs <- datadea$ud_outputs
    nuo <- length(ud_outputs)
    
    if (!all(ud_outputs %in% 1:no)) {
      stop("Invalid set of undesirable outputs.")
    }
    
    if (is.null(vtrans_o)) {
      vtrans_o <- rep(NA, nuo)
    } else if ((length(vtrans_o) == 1) && (nuo > 1)) {
      vtrans_o <- rep(vtrans_o, nuo)
    } else if (length(vtrans_o) != nuo) {
      stop("Translation vector vtrans_o must be NULL, a constant or a vector of the
           same length as ud_outputs.")
    }
    
    if (nuo > 0) {
      for (i in 1:nuo) {
        if (is.na(vtrans_o[i])) {
          vtrans_o[i] <- max(output$mR[ud_outputs[i], ] + output$dR[ud_outputs[i], ]) + 1
        }
      }
    }
    
    i <- 1
    for (j in ud_outputs) {
      
      aux <- output$dL[j, ]
      output$dL[j, ] <- output$dR[j, ]
      output$dR[j, ] <- aux
      
      aux <- output$mR[j, ]
      output$mR[j, ] <- vtrans_o[i] - output$mL[j, ]
      output$mL[j, ] <- vtrans_o[i] - aux
      
      if (length(output$mL[j, (output$mL[j, ] - output$dL[j, ]) < 0] > 0)) {
        stop("There is at least one negative bad output. Change translation vector vtrans_o.")
        print(j)
      }
      
      i <- i + 1
    }
    
    u_datadea <- datadea
    u_datadea$input <- input
    u_datadea$output <- output
    
    return(list(
      u_datadea = structure(u_datadea, class = "deadata_fuzzy"),
      vtrans_i = vtrans_i,
      vtrans_o = vtrans_o))
    
  } else {
    stop("Data should be of class deadata or deadata_fuzzy.
         Run read_data or read_data_fuzzy function first!")
  }
    
}