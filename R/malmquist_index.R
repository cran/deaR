#' @title Malmquist index
#'   
#' @description This function calculates the input/output oriented
#' Malmquist productivity index under constant or variable returns-to-scale.
#' 
#' @note In the results: EC = Efficiency Change, PTEC = Pure Technical Efficiency Change,
#' SEC = Scale Efficiency Change, TC = Technological Change, MI = Malmquist Index 
#' 
#' @usage malmquist_index(datadealist,
#'                 dmu_eval = NULL,
#'                 dmu_ref = NULL,
#'                 orientation = c("io", "oo"),
#'                 rts = c("crs", "vrs"),
#'                 type1 = c("cont", "seq", "glob"),
#'                 type2 = c("fgnz", "rd", "gl", "bias"),
#'                 tc_vrs = FALSE,
#'                 vtrans_i = NULL,
#'                 vtrans_o = NULL)
#' 
#' @param datadealist A list with the data at different times, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param orientation A string, equal to "io" (input oriented) or "oo" (output oriented).
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant) or
#'            "vrs" (variable).
#' @param type1 A string, equal to "cont" (contemporary), "seq" (sequential) or "glob"
#' (global).
#' @param type2 A string, equal to "fgnz" (Fare et al. 1994), "rd" (Ray and Desli 1997),
#' "gl" (generalized) or "bias" (biased).
#' @param tc_vrs Logical. If it is \code{FALSE}, it computes the vrs bias malmquist index by using
#' the technical change under crs (Fare and Grosskopf 1996). Otherwise, it uses the technical
#' change under vrs.
#' @param vtrans_i Numeric vector of translation for undesirable inputs in non-directional
#' basic models. If \code{vtrans_i[i]} is \code{NA}, then it applies the "max + 1" translation
#' to the i-th undesirable input. If \code{vtrans_i} is a constant, then it applies
#' the same translation to all undesirable inputs. If \code{vtrans_i} is \code{NULL},
#' then it applies the "max + 1" translation to all undesirable inputs.
#' @param vtrans_o Numeric vector of translation for undesirable outputs in
#' non-directional basic models, analogous to \code{vtrans_i}, but applied to outputs.
#'   
#' @return A numeric list with Malmquist index and other parameters.
#' 
#' @author 
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#' 
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'  
#' @references 
#' Caves, D.W.; Christensen, L. R.; Diewert, W.E. (1982). “The Economic Theory of
#' Index Numbers and the Measurement of Input, Output, and Productivity”. Econometrica,
#' 50(6), 1393-1414. 
#'  
#' Fare, R.; Grifell-Tatje, E.; Grosskopf, S.; Lovell, C.A.K. (1997). "Biased Technical
#' Change and the Malmquist Productivity Index". Scandinavian Journal of Economics,
#' 99(1), 119-127.
#' 
#' Fare, R.; Grosskopf, S.; Lindgren, B.; Roos, P. (1989). “Productivity Developments
#' in Swedish Hospitals: A Malmquist Output Index Approach”. Discussion paper n. 89-3.
#' Southern Illinois University. Illinois.
#' 
#' Fare, R.; Grosskopf, S.; Lindgren, B.; Roos, P. (1992). “Productivity changes
#' in Swedish Pharmacies 1980-89: A nonparametric Malmquist Approach”. Journal of
#' productivity Analysis, 3(3), 85-101. 
#' 
#' Fare, R.; Grosskopf, S.; Norris, M.; Zhang, Z. (1994). “Productivity Growth,
#' Technical Progress, and Efficiency Change in Industrialized Countries”.
#' American Economic Review, 84(1), 66-83. 
#'
#' Fare, R.; Grosskopf, S.; Roos, P. (1998), Malmquist Productivity Indexes: A Survey
#' of Theory and Practice. In: Fare R., Grosskopf S., Russell R.R. (eds)
#' Index Numbers: Essays in Honour of Sten Malmquist. Springer.
#' 
#' Grifell-Tatje, E.; Lovell, C.A.K. (1999). "A Generalized Malmquist productivity index".
#' Top, 7(1), 81-101.  
#'
#' Pastor, J.T.; Lovell, C.A.k. (2005). "A global Malmquist productiviyt index".
#' Economics Letters, 88, 266-271.  
#' 
#' Ray, S.C.; Desli, E. (1997). "Productivity Growth, Technical Progress, and
#' Efficiency Change in Industrialized Countries: Comment". The American Economic Review,
#' 87(5), 1033-1039.
#'
#' Shestalova, V. (2003). "Sequential Malmquist Indices of Productivity Growth:
#' An Application to OECD Industrial Activities". Journal of Productivity Analysis,
#' 19, 211-226.
#'
#'
#' @examples 
#' # Example 1. With dataset in wide format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("Economy")
#' data_example <- read_malmquist(datadea = Economy,
#'                                nper = 5, 
#'                                arrangement = "horizontal",
#'                                ni = 2, 
#'                                no = 1)
#' result <- malmquist_index(data_example, orientation = "io")
#' mi <- result$mi
#' effch <- result$ec
#' tech <- result$tc
#' 
#' # Example 2. With dataset in long format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("EconomyLong")
#' data_example2 <- read_malmquist(EconomyLong,
#'                                 percol = 2, 
#'                                 arrangement = "vertical",
#'                                 inputs = 3:4, 
#'                                 outputs = 5)
#' result2 <- malmquist_index(data_example2, orientation = "io")
#' mi2 <- result2$mi
#' effch2 <- result2$ec
#' tech2 <- result2$tc
#' 
#' # Example 3. Replication of results in Grifell-Tatje and Lovell (1999, p. 100).
#' data("Grifell_Lovell_1999")
#' data_example <- read_malmquist(Grifell_Lovell_1999,
#'                                percol = 1,
#'                                dmus = 2,
#'                                inputs = 3,
#'                                outputs = 4,
#'                                arrangement = "vertical")
#' result_fgnz <- malmquist_index(data_example,
#'                                orientation = "oo",
#'                                rts = "vrs",
#'                                type1 = "cont",
#'                                type2 = "fgnz")
#' mi_fgnz <- result_fgnz$mi 
#' 
#' result_rd <- malmquist_index(data_example,
#'                              orientation = "oo",
#'                              rts = "vrs",
#'                              type1 = "cont",
#'                              type2 = "rd")
#' mi_rd <- result_rd$mi
#'  
#' result_gl <- malmquist_index(data_example,
#'                              orientation = "oo",
#'                              rts = "vrs",
#'                              type1 = "cont",
#'                              type2 = "gl")
#' mi_gl <- result_gl$mi
#'                               
#' @import lpSolve
#' 
#' @export
  
malmquist_index <- function(datadealist,
                            dmu_eval = NULL,
                            dmu_ref = NULL,
                            orientation = c("io", "oo"),
                            rts = c("crs", "vrs"),
                            type1 = c("cont", "seq", "glob"),
                            type2 = c("fgnz", "rd", "gl", "bias"),
                            tc_vrs = FALSE,
                            vtrans_i = NULL,
                            vtrans_o = NULL) {
  
  nt <- length(datadealist)
  
  if (nt < 2) {
    stop("Number of times should be >= 2.")
  }
  
  # Cheking data
  for (t in 1:nt) {
    if (!is.deadata(datadealist[[t]])) {
      stop("Data should be of class deadata. Run read_data function first!")
    }
  }
  
  pernames <- names(datadealist)
  minames <- pernames[-1]
  
  dmunames <- datadealist[[1]]$dmunames
  nd <- length(dmunames) # number of dmus
  for (t in 2:nt) {
    if (nd != length(datadealist[[t]]$dmunames)) {
      stop("Number of DMUs does not coincide.")
    }
  }
  
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
  
  ni <- nrow(datadealist[[1]]$input) # number of inputs
  no <- nrow(datadealist[[1]]$output) # number of outputs
  for (t in 2:nt) {
    if (ni != nrow(datadealist[[t]]$input)) {
      stop("Number of inputs does not coincide.")
    }
    if (no != nrow(datadealist[[t]]$output)) {
      stop("Number of outputs does not coincide.")
    }
  }
  
  # Checking types
  type1 <- tolower(type1)
  type1 <- match.arg(type1)
  type2 <- tolower(type2)
  type2 <- match.arg(type2)
  
  # Checking rts
  rts <- tolower(rts)
  rts <- match.arg(rts)
  
  # Checking orientation
  orientation <- tolower(orientation)
  orientation <- match.arg(orientation)
  
  # Checking undesirable variables
  if (!is.null(datadealist[[1]]$ud_inputs) || !is.null(datadealist[[1]]$ud_outputs)) {
    warning("There are undesirable variables and some efficiencies are computed under constant
            returns to scale.")
  }

  input <- array(0, dim = c(ni, nd, nt))
  output <- array(0, dim = c(no, nd, nt))
  for (t in 1:nt) {
    input[, , t] <- datadealist[[t]]$input
    output[, , t] <- datadealist[[t]]$output
  }
  
  mi <- matrix(0, nrow = nt - 1, ncol = nde)
  colnames(mi) <- dmunames[dmu_eval]
  rownames(mi) <- minames
  eff <- matrix(0, nrow = nt, ncol = nde) # efficiencies crs
  colnames(eff) <- dmunames[dmu_eval]
  rownames(eff) <- pernames
  effv <- eff # efficiencies vrs
  
  if (type1 == "glob") {
    type2 <- NULL
    effg <- eff # efficiencies with global frontier crs 
    effgv <- eff # efficiencies with global frontier vrs 
    if (rts == "vrs") {
      eff12 <- mi # DMU forward
      effv12 <- mi # DMU forward vrs
    }
  } else {
    if (type2 == "fgnz") {
      eff12 <- mi # DMU forward
      eff21 <- mi # Frontier forward
    } else if (type2 == "rd") {
      eff12 <- mi # DMU forward
      effv12 <- mi # DMU forward vrs
      effv21 <- mi # Frontier forward vrs
    } else if (type2 == "gl") {
      eff12y <- mi # DMU input forward
      effv12 <- mi # DMU forward vrs
      effv12y <- mi # DMU input forward vrs
    } else if (type2 == "bias") {
      eff12 <- mi # DMU forward
      eff21 <- mi # Frontier forward
      eff12y <- mi # DMU input forward
      eff22y <- mi # Frontier and DMU input forward
      effv12 <- mi # DMU forward
      effv21 <- mi # Frontier forward
      effv12y <- mi # DMU input forward
      effv22y <- mi # Frontier and DMU input forward
    }
  }
  
  if (type1 == "cont") {
    
    datadea <- datadealist[[1]]
    eff[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    effv[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, rts = "vrs",
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    
    dmu_ref_cont <- dmu_ref + 1
    
    for (t in 2:nt) {
      
      datadea <- datadealist[[t]]
      eff[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      effv[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, rts = "vrs",
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      
      # Intertemporal scores 
      
      for (i in 1:nde) {
        
        ii <- dmu_eval[i]
        datadea$dmunames <- c(dmunames[ii], dmunames)
        
        if (type2 == "fgnz") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t - 1], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        matrix(input[, , t], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         matrix(output[, , t], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff21[t - 1, i] <- 
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "rd") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t - 1], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        matrix(input[, , t], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         matrix(output[, , t], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          effv21[t - 1, i] <- 
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "gl") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t - 1], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "bias") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t - 1], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        matrix(input[, , t], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         matrix(output[, , t], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          eff22y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv22y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        }
        
      }
    }

  } else if (type1 == "seq") {
    
    datadea <- datadealist[[1]]
    eff[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    effv[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, rts = "vrs",
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    
    inputfront1 <- matrix(input[, , 1], nrow = ni, ncol = nd)
    outputfront1 <- matrix(output[, , 1], nrow = no, ncol = nd)
    dmu_ref_seq1 <- dmu_ref + 1
    
    for (t in 2:nt) {
      
      inputfront2 <- cbind(inputfront1, matrix(input[, , t], nrow = ni, ncol = nd))
      outputfront2 <- cbind(outputfront1, matrix(output[, , t], nrow = no, ncol = nd))
      dmu_ref_seq2 <- c(dmu_ref_seq1, dmu_ref + (t - 1) * nd + 1)
      dmu_eval_seq <- dmu_eval + (t - 1) * nd
      ncol1 <- (t - 1) * nd + 1
      ncol2 <- t * nd + 1
      
      datadea <- datadealist[[t]]
      
      datadea$dmunames <- rep(dmunames, t)
      datadea$input <- inputfront2
      datadea$output <- outputfront2
      eff[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval_seq, dmu_ref = dmu_ref_seq2 - 1,
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      effv[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval_seq, dmu_ref = dmu_ref_seq2 - 1, rts = "vrs",
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      
      # Intertemporal scores 
      
      for (i in 1:nde) {
        
        ii <- dmu_eval[i]
        dmunames_seq1 <- c(dmunames[ii], rep(dmunames, t - 1))
        dmunames_seq2 <- c(dmunames_seq1, dmunames)
        
        if (type2 == "fgnz") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        inputfront1),
                                  nrow = ni, ncol = ncol1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          datadea$dmunames <- dmunames_seq1
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        inputfront2),
                                  nrow = ni, ncol = ncol2)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         outputfront2),
                                   nrow = no, ncol = ncol2)
          datadea$dmunames <- dmunames_seq2
          eff21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "rd") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        inputfront1),
                                  nrow = ni, ncol = ncol1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          datadea$dmunames <- dmunames_seq1
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        inputfront2),
                                  nrow = ni, ncol = ncol2)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         outputfront2),
                                   nrow = no, ncol = ncol2)
          datadea$dmunames <- dmunames_seq2
          effv21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "gl") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        inputfront1),
                                  nrow = ni, ncol = ncol1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          datadea$dmunames <- dmunames_seq1
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          eff12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        } else if (type2 == "bias") {
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        inputfront1),
                                  nrow = ni, ncol = ncol1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          datadea$dmunames <- dmunames_seq1
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         outputfront1),
                                   nrow = no, ncol = ncol1)
          eff12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq1, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t - 1], nrow = ni, ncol = 1),
                                        inputfront2),
                                  nrow = ni, ncol = ncol2)
          datadea$output <- matrix(cbind(matrix(output[, ii, t - 1], nrow = no, ncol = 1),
                                         outputfront2),
                                   nrow = no, ncol = ncol2)
          datadea$dmunames <- dmunames_seq2
          eff21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv21[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        inputfront2),
                                  nrow = ni, ncol = ncol2)
          eff22y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv22y[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_seq2, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        }
        
      }
      
      inputfront1 <- inputfront2
      outputfront1 <- outputfront2
      dmu_ref_seq1 <- dmu_ref_seq2
    }
    
  } else if (type1 == "glob") {
    
    # Frontier usual
    
    datadea <- datadealist[[1]]
    eff[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    effv[1, ] <- efficiencies(
      model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, rts = "vrs",
                  orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                  vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
    )
    
    dmu_ref_cont <- dmu_ref + 1
    
    for (t in 2:nt) {
      
      datadea <- datadealist[[t]]
      eff[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref,
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      effv[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval, dmu_ref = dmu_ref, rts = "vrs",
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      
      # Intertemporal scores 
      
      if (rts == "vrs") {
        
        for (i in 1:nde) {
          
          ii <- dmu_eval[i]
          datadea$dmunames <- c(dmunames[ii], dmunames)
          
          datadea$input <- matrix(cbind(matrix(input[, ii, t], nrow = ni, ncol = 1),
                                        matrix(input[, , t - 1], nrow = ni, ncol = nd)),
                                  nrow = ni, ncol = nd + 1)
          datadea$output <- matrix(cbind(matrix(output[, ii, t], nrow = no, ncol = 1),
                                         matrix(output[, , t - 1], nrow = no, ncol = nd)),
                                   nrow = no, ncol = nd + 1)
          eff12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont,
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
          effv12[t - 1, i] <-
            model_basic(datadea, dmu_eval = 1, dmu_ref = dmu_ref_cont, rts = "vrs",
                        orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                        vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)$DMU[[1]]$efficiency
          
        }
        
      }
      
    }
    
    # Frontier global
    
    datadea <- datadealist[[1]]
    inputfront <- matrix(0, nrow = ni, ncol = nt * nd)
    outputfront <- matrix(0, nrow = no, ncol = nt * nd)
    dmu_ref_glob <- rep(0, nt * ndr)
    for (t in 1:nt) {
      inputfront[, ((t - 1) * nd + 1):(t * nd)] <- matrix(input[, , t], nrow = ni)
      outputfront[, ((t - 1) * nd + 1):(t * nd)] <- matrix(output[, , t], nrow = no)
      dmu_ref_glob[((t - 1) * ndr + 1):(t * ndr)] <- dmu_ref + (t - 1) * nd
    }
    datadea$input <- inputfront
    datadea$output <- outputfront
    datadea$dmunames <- rep(dmunames, nt)
    
    for (t in 1:nt) {
      
      dmu_eval_glob <- dmu_eval + (t - 1) * nd
      effg[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval_glob, dmu_ref = dmu_ref_glob,
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      effgv[t, ] <- efficiencies(
        model_basic(datadea, dmu_eval = dmu_eval_glob, dmu_ref = dmu_ref_glob, rts = "vrs",
                    orientation = orientation, maxslack = FALSE, compute_target = FALSE,
                    vtrans_i = vtrans_i, vtrans_o = vtrans_o, silent_ud = TRUE)
      )
      
    }
    
  }
  
  if (orientation == "oo") {
    eff <- 1 / eff
    effv <- 1 / effv
    if (type1 == "glob") {
      effg <- 1 / effg 
      effgv <- 1 / effgv  
      if (rts == "vrs") {
        eff12 <- 1 / eff12
        effv12 <- 1 / effv12
      }
    } else {
      if (type2 == "fgnz") {
        eff12 <- 1 / eff12
        eff21 <- 1 / eff21
      } else if (type2 == "rd") {
        eff12 <- 1 / eff12
        effv12 <- 1 / effv12
        effv21 <- 1 / effv21
      } else if (type2 == "gl") {
        eff12y <- 1 / eff12y
        effv12 <- 1 / effv12
        effv12y <- 1 / effv12y
      } else if (type2 == "bias") {
        eff12 <- 1 / eff12
        eff21 <- 1 / eff21
        eff12y <- 1 / eff12y 
        eff22y <- 1 / eff22y
        effv12 <- 1 / effv12
        effv21 <- 1 / effv21
        effv12y <- 1 / effv12y 
        effv22y <- 1 / effv22y
      }
    }
  }
  
  ec <- NULL
  pech <- NULL
  sech <- NULL
  obtech <- NULL
  ibtech <- NULL
  matech <- NULL
  
  if (type1 == "glob") {
    if (rts == "crs") {
      ec <- eff[-1, ] / eff[-nt, ]
      tc <- (effg[-1, ] * eff[-nt, ]) / (eff[-1, ] * effg[-nt, ])
      mi <- ec * tc
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency.glob.crs = effg, efficiency.glob.vrs = effgv)
    } else {
      pech <- effv[-1, ] / effv[-nt, ]
      tc <- (effgv[-1, ] * effv[-nt, ]) / (effv[-1, ] * effgv[-nt, ])
      sech <- (effv[-nt, ] / eff[-nt, ]) / (effv12 / eff12) # Grifell-Tatjé and Lovell (1999)
      mi <- tc * pech * sech
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency.glob.crs = effg, efficiency.glob.vrs = effgv,
                      efficiency_t_t1.crs = eff12, efficiency_t_t1.vrs = effv12)
    }
  } else {
    
    if (type2 == "fgnz") {
      # mi <- sqrt((eff12 * eff[-1, ]) / (eff[-nt, ] * eff21))
      tc <- ((eff12 / eff[-1, ]) * (eff[-nt, ] / eff21)) ^ 0.5
      if (rts == "crs") {
        ec <- eff[-1, ] / eff[-nt, ]
        mi <- ec * tc
      } else {
        pech <- effv[-1, ] / effv[-nt, ]
        sech <- (effv[-nt, ] / eff[-nt, ]) / (effv[-1, ] / eff[-1, ])
        mi <- tc * pech * sech
      }
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency_t_t1.crs = eff12, efficiency_t1_t.crs = eff21)
    } else if (type2 == "rd") {
      tc <- effv12 / effv[-1, ] # Grifell-Tatjé and Lovell (1999)
      # tc <- ((effv[-nt, ] / effv21) * (effv12 / effv[-1, ])) ^ 0.5  # Ray and Desli (1997)
      if (rts == "crs") {
        warning("Descomposition only under variable returns-to-scale (vrs). Parameter rts has
                been set to vrs.")
        rts <- "vrs"
      }
      pech <- effv[-1, ] / effv[-nt, ]
      sech <- (effv[-nt, ] / eff[-nt, ]) / (effv12 / eff12) # Grifell-Tatjé and Lovell (1999)
      # sech <- (((effv[-nt, ] / eff[-nt, ]) / (effv12 / eff12)) * ((effv21 / eff21) / (effv[-1, ] / (eff[-1, ])))) ^ 0.5   # Ray and Desli (1997)
      mi <- tc * pech * sech
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency_t_t1.crs = eff12, efficiency_t_t1.vrs = effv12,
                      efficiency_t1_t.vrs = effv21)
    } else if (type2 == "gl") { # generalized
      tc <- effv12 / effv[-1, ]
      if (rts == "crs") {
        warning("Descomposition only under variable returns-to-scale (vrs). Parameter rts has
                been set to vrs.")
        rts <- "vrs"
      }
      pech <- effv[-1, ] / effv[-nt, ]
      sech <- (effv[-nt, ] / eff[-nt, ]) / (effv12y / eff12y) 
      mi <- tc * pech * sech
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency_t_xt1.crs = eff12y, efficiency_t_t1.vrs = effv12,
                      efficiency_t_xt1.vrs = effv12y)
    } else if (type2 == "bias") {
      if (rts == "crs") {
        obtech <- sqrt((eff12 * eff22y) / (eff[-1, ] * eff12y))
        ibtech <- sqrt((eff21 * eff12y) / (eff[-nt, ] * eff22y))
        matech <- eff[-nt, ] / eff21
        tc <- obtech * ibtech * matech
        ec <- eff[-1, ] / eff[-nt, ]
        mi <- ec * tc
      } else {
        if (tc_vrs) {
          obtech <- sqrt((effv12 * effv22y) / (effv[-1, ] * effv12y))
          ibtech <- sqrt((effv21 * effv12y) / (effv[-nt, ] * effv22y))
          matech <- effv[-nt, ] / effv21
          tc <- obtech * ibtech * matech
        } else {
          warning("By default, technical change (tc) is measured relative to constant returns to scale (crs)
                   (Fare and Grosskopf 1996). For tc under variable returns to scale (vrs)
                   set the logical parameter tc_vrs to TRUE.")
          obtech <- sqrt((eff12 * eff22y) / (eff[-1, ] * eff12y))
          ibtech <- sqrt((eff21 * eff12y) / (eff[-nt, ] * eff22y))
          matech <- eff[-nt, ] / eff21
          tc <- obtech * ibtech * matech
        }
        pech <- effv[-1, ] / effv[-nt, ]
        sech <- (eff[-1, ] / effv[-1, ]) / (eff[-nt, ] / effv[-nt, ]) #fgnz
        mi <- pech * sech * tc
      }
      eff_all <- list(efficiency.crs = eff, efficiency.vrs = effv,
                      efficiency_t_t1.crs = eff12, efficiency_t1_t.crs = eff21,
                      efficiency_t_xt1.crs = eff12y, efficiency_t1_xt1.crs = eff22y,
                      efficiency_t_t1.vrs = effv12, efficiency_t1_t.vrs = effv21,
                      efficiency_t_xt1.vrs = effv12y, efficiency_t1_xt1.vrs = effv22y)
    }
    
  }
  
  minames <- pernames[-1]
  rownames(mi) <- minames
  rownames(tc) <- minames

  if (!is.null(ec)) {
    if (!is.matrix(ec)) {
      ec <- matrix(ec, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(ec) <- minames
    }
  }
  if (!is.null(pech)) {
    if (!is.matrix(pech)) {
      pech <- matrix(pech, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(pech) <- minames
    }
  }
  if (!is.null(sech)) {
    if (!is.matrix(sech)) {
      sech <- matrix(sech, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(sech) <- minames
    }
  }
  if (!is.null(obtech)) {
    if (!is.matrix(obtech)) {
      obtech <- matrix(obtech, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(obtech) <- minames
    }
  }
  if (!is.null(ibtech)) {
    if (!is.matrix(ibtech)) {
      ibtech <- matrix(ibtech, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(ibtech) <- minames
    }
  }
  if (!is.null(matech)) {
    if (!is.matrix(matech)) {
      matech <- matrix(matech, nrow = nt - 1, ncol = nde, dimnames = list(minames, names(dmu_eval)))
    } else {
      rownames(matech) <- minames
    }
  }
  
  deaOutput <- list(mi = mi,
                    ec = ec,
                    tc = tc,
                    pech = pech,
                    sech = sech,
                    obtech = obtech,
                    ibtech = ibtech,
                    matech = matech,
                    eff_all = eff_all,
                    datadealist = datadealist,
                    dmu_eval = dmu_eval,
                    dmu_ref = dmu_ref,
                    orientation = orientation,
                    rts = rts,
                    type1 = type1,
                    type2 = type2,
                    tc_vrs = tc_vrs,
                    vtrans_i = vtrans_i,
                    vtrans_o = vtrans_o, silent_ud = TRUE,
                    modelname = "malmquist")
  return(structure(deaOutput, class = "dea"))
  
}