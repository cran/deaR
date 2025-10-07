#' @title read_data
#'  
#' @description This function is deprecated. Use \code{make_deadata} instead.
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
#' @param dmus Column (number or name) of DMUs (optional). By default, it is the
#' first column. If there is not any DMU column, then it must be \code{NULL}.
#' @param ni Number of inputs, if inputs are in columns 2:(\code{ni} + 1) (if DMUs
#' are in the first column) or 1:\code{ni} (no DMUs column).
#' @param no Number of outputs, if outputs are in columns (\code{ni} + 2):(\code{ni} +
#' \code{no} + 1) (if DMUs are in the first column) or (\code{ni} + 1):(\code{ni} +
#' \code{no}) (no DMUs column). If not specified, DMUs are in the first column.
#' @param inputs Columns (numbers or names) of inputs (optional). It prevails over \code{ni}.
#' Alternatively to \code{datadea}, it can be a matrix with the inputs (DMUs in columns).
#' In this case, DMUs names are taken from the columns names.
#' @param outputs Columns (numbers or names) of outputs (optional). It prevails over \code{no}.
#' Alternatively to \code{datadea}, it can be a matrix with the outputs (DMUs in columns).
#' @param nc_inputs A numeric vector containing the indices of non-controllable inputs.
#' @param nc_outputs A numeric vector containing the indices of non-controllable outputs.
#' @param nd_inputs A numeric vector containing the indices of non-discretionary inputs.
#' @param nd_outputs A numeric vector containing the indices of non-discretionary outputs.
#' @param ud_inputs A numeric vector containing the indices of undesirable (good) inputs.
#' @param ud_outputs A numeric vector containing the indices of undesirable (bad) outputs.
#' 
#' @returns An object of class \code{deadata}. Any function \code{model_*} can be applied.
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
  
  .Deprecated("make_deadata")
  
  make_deadata(datadea = datadea,
               ni = ni,
               no = no,
               dmus = dmus,
               inputs = inputs,
               outputs = outputs,
               nc_inputs = nc_inputs,
               nc_outputs = nc_outputs,
               nd_inputs = nd_inputs,
               nd_outputs = nd_outputs,
               ud_inputs = ud_inputs,
               ud_outputs = ud_outputs)
  
}