#' @title read_data_fuzzy
#'  
#' @description This function is deprecated. Use \code{make_deadata_fuzzy} instead.
#' 
#' @usage read_data_fuzzy(datadea,
#'                 dmus = 1,
#'                 inputs.mL = NULL,
#'                 inputs.mR = NULL,
#'                 inputs.dL = NULL,
#'                 inputs.dR = NULL,
#'                 outputs.mL = NULL,
#'                 outputs.mR = NULL,
#'                 outputs.dL = NULL,
#'                 outputs.dR = NULL,
#'                 nc_inputs = NULL,
#'                 nc_outputs = NULL,
#'                 nd_inputs = NULL,
#'                 nd_outputs = NULL,
#'                 ud_inputs = NULL,
#'                 ud_outputs = NULL)
#'          
#' @param datadea Data frame with DEA data.
#' @param dmus Column (number or name) of DMUs (optional). By default, it is the first
#' column. If there is not any DMU column, then it must be \code{NULL}.
#' @param inputs.mL Where are (columns) the \code{mL} (left centers) of trapezoidal fuzzy
#' inputs in \code{datadea}. If an input is triangular or crisp, we put the column
#' where the centers or the crisp values are, respectively.
#' 
#' Alternatively to \code{datadea}, \code{inputs.mL} can be a matrix of size (number of inputs x
#' number of DMUs) with the \code{mL} of trapezoidal fuzzy inputs, the centers of
#' triangular inputs, and the crisp values of crisp inputs. In this case, DMUs names are
#' taken from the columns names.
#' @param inputs.mR Where are (columns) the \code{mR} (right centers) of trapezoidal fuzzy
#' inputs in \code{datadea}. If an input is triangular or crisp, we put \code{NA}.
#' 
#' Alternatively to \code{datadea}, \code{inputs.mR} can be a matrix of size (number of inputs x
#' number of DMUs) with the \code{mR} of trapezoidal fuzzy inputs, the centers of
#' triangular inputs, and the crisp values of crisp inputs. If all inputs are triangular or
#' crisp, then \code{inputs.mR} must be NULL (default) or equal to \code{inputs.mL}.
#' @param inputs.dL Where are (columns) the \code{dL} (left radii) of trapezoidal and
#' triangular fuzzy inputs in \code{datadea}. If an input is symmetric, we put the column
#' where the radii are. If an input is rectangular or crisp, we put \code{NA}.
#' 
#' Alternatively to \code{datadea}, \code{inputs.dL} can be a matrix of size (number of inputs x
#' number of DMUs) with the \code{dL} of trapezoidal and triangular fuzzy inputs. If an
#' input is rectangular or crisp, its radius is zero. If all inputs are rectangular or
#' crisp, then \code{inputs.dL} must be NULL (default) or a zero matrix.
#' @param inputs.dR Where are (columns) the \code{dR} (right radii) of trapezoidal and
#' triangular fuzzy inputs in \code{datadea}. If an input is symmetric, rectangular or
#' crisp, we put \code{NA}.
#' 
#' Alternatively to \code{datadea}, \code{inputs.dR} can be a matrix of size (number of inputs x
#' number of DMUs) with the \code{dR} of trapezoidal and triangular fuzzy inputs. If an
#' input is rectangular or crisp, its radius is zero. If all inputs are symmetric,
#' rectangular or crisp, then \code{inputs.dR} must be NULL (default) or equal to
#' \code{inputs.dL}.
#' @param outputs.mL Analogous to \code{inputs.mL}, but relating to outputs.
#' @param outputs.mR Analogous to \code{inputs.mR}, but relating to outputs.
#' @param outputs.dL Analogous to \code{inputs.dL}, but relating to outputs.
#' @param outputs.dR Analogous to \code{inputs.dR}, but relating to outputs.
#' @param nc_inputs A numeric vector containing the indices of non-controllable inputs.
#' @param nc_outputs A numeric vector containing the indices of non-controllable outputs.
#' @param nd_inputs A numeric vector containing the indices of non-discretionary inputs.
#' @param nd_outputs A numeric vector containing the indices of non-discretionary outputs.
#' @param ud_inputs A numeric vector containing the indices of undesirable (good) inputs.
#' @param ud_outputs A numeric vector containing the indices of undesirable (bad) outputs.
#' 
#' @returns An object of class \code{deadata_fuzzy}. Any function \code{modelfuzzy_*}
#' can be applied.
#'                                 
#' @export

read_data_fuzzy <- function(datadea,
                            dmus = 1,
                            inputs.mL = NULL,
                            inputs.mR = NULL,
                            inputs.dL = NULL,
                            inputs.dR = NULL,
                            outputs.mL = NULL,
                            outputs.mR = NULL,
                            outputs.dL = NULL,
                            outputs.dR = NULL,
                            nc_inputs = NULL,
                            nc_outputs = NULL,
                            nd_inputs = NULL,
                            nd_outputs = NULL,
                            ud_inputs = NULL,
                            ud_outputs = NULL) {
  
  .Deprecated("make_deadata_fuzzy")
  
  make_deadata_fuzzy(datadea = datadea,
                     dmus = dmus,
                     inputs.mL = inputs.mL,
                     inputs.mR = inputs.mR,
                     inputs.dL = inputs.dL,
                     inputs.dR = inputs.dR,
                     outputs.mL = outputs.mL,
                     outputs.mR = outputs.mR,
                     outputs.dL = outputs.dL,
                     outputs.dR = outputs.dR,
                     nc_inputs = nc_inputs,
                     nc_outputs = nc_outputs,
                     nd_inputs = nd_inputs,
                     nd_outputs = nd_outputs,
                     ud_inputs = ud_inputs,
                     ud_outputs = ud_outputs)
  
}