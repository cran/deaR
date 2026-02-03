#' @title Fuzzy DEA model.
#'   
#' @description Solve the fuzzy DEA model by Kao and Liu (2000)
#' 
#' @usage modelfuzzy_kaoliu(datadea,
#'                   dmu_eval = NULL,
#'                   kaoliu_modelname = c("basic", "additive", "addmin",
#'                                        "addsupereff", "deaps", "dir", "fdh",
#'                                        "lgo", "multiplier", "nonradial",
#'                                        "profit", "qgo", "rdm", "sbmeff",
#'                                        "sbmsupereff", "supereff"),
#'                   alpha = 1,
#'                   ...)
#' 
#' @param datadea A \code{deadata_fuzzy} object, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param kaoliu_modelname a string containing the name of the model.
#' @param alpha A numeric vector with the alpha-cuts (in [0,1]). If \code{alpha}>1,
#' it determines the number of alpha-cuts, equispatially distributed in [0,1]. 
#' @param ... \code{dmu_ref}, \code{orientation}, \code{rts} and other model parameters.
#'   
#' @returns A list of class \code{dea_fuzzy} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{orientation}, \code{rts}, \code{alpha},
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
#' Boscá, J.E.; Liern, V.; Sala, R.; Martínez, A. (2011). "Ranking Decision Making
#' Units by Means of Soft Computing DEA Models". International Journal of Uncertainty,
#' Fuzziness and Knowledge-Based Systems, 19(1), p.115-134. 
#' 
#' Emrouznejad, A.; Tavana, M.; Hatami-Marbini, A. (2014). “The State of the Art
#' in Fuzzy Data Envelopment Analysis”, in A. Emrouznejad and M. Tavana (eds.),
#' Performance Measurement with Fuzzy Data Envelopment Analysis. Studies in Fuzziness
#' and Soft Computing 309. Springer, Berlin.
#' 
#' Hatami-Marbini, A.; Emrouznejad, A.; Tavana, M. (2011). "A Taxonomy and Review
#' of the Fuzzy Data Envelopment Analysis Literature: Two Decades in the Making",
#' European Journal of Operational Research, 214, 457–472.
#' 
#' Kao, C.; Liu, S.T. (2000). “Fuzzy efficiency measures in data envelopment analysis,
#' Fuzzy Sets and Systems”, 119, 149–160.
#' 
#' Kao, C., Liu, S.T., (2000). “Data envelopment analysis with missing data:
#' An application to university libraries in Taiwan”, Journal of the Operational
#' Research Society, 51, 897–905. 
#' 
#' Kao, C., Liu, S.T. (2003). “A mathematical programming approach to fuzzy
#' efficiency ranking”, International Journal of Production Economics, 85.
#' 
#' @examples 
#' # Example 1. 
#' # Replication of results in Boscá, Liern, Sala and Martínez (2011, p.125)
#' data("Leon2003")
#' data_example <- make_deadata_fuzzy(datadea = Leon2003,
#'                                    inputs.mL = 2, 
#'                                    inputs.dL = 3, 
#'                                    outputs.mL = 4, 
#'                                    outputs.dL = 5)
#' result <- modelfuzzy_kaoliu(data_example,
#'                             kaoliu_modelname = "basic", 
#'                             alpha = seq(0, 1, by = 0.1), 
#'                             orientation = "io", 
#'                             rts = "vrs")
#' efficiencies(result)
#' 
#' # Example 2.
#' # Replication of results in Kao and Liu (2003, p.152)
#' data("Kao_Liu_2003")
#' data_example <- make_deadata_fuzzy(Kao_Liu_2003, 
#'                                    inputs.mL = 2, 
#'                                    outputs.mL = 3:7, 
#'                                    outputs.dL = c(NA, NA, 8, NA, 10),
#'                                    outputs.dR = c(NA, NA, 9, NA, 11))
#' result <- modelfuzzy_kaoliu(data_example, 
#'                             kaoliu_modelname = "basic", 
#'                             orientation = "oo", 
#'                             rts = "vrs", 
#'                             alpha = 0)
#' sol <- efficiencies(result)
#' eff <- data.frame(1 / sol$Worst, 1 / sol$Best)
#' names(eff) <- c("eff_lower", "eff_upper")
#' eff
#' 
#' @seealso \code{\link{model_basic}}, \code{\link{model_multiplier}},
#' \code{\link{modelfuzzy_possibilistic}}, \code{\link{modelfuzzy_guotanaka}}
#' 
#' @import lpSolve
#' 
#' @export
  
modelfuzzy_kaoliu <-
  function(datadea,
           dmu_eval = NULL,
           kaoliu_modelname = c("basic", "additive", "addmin", "addsupereff", "deaps",
                                "dir", "fdh", "lgo", "multiplier", "nonradial", "profit",
                                "qgo", "rdm", "sbmeff", "sbmsupereff", "supereff"),
           alpha = 1,
           ...) {
 
  # Checking whether datadea is of class "deadata_fuzzy" or not...  
  if (!is.deadata_fuzzy(datadea)) {
    stop("Data should be of class deadata_fuzzy. Run make_deadata_fuzzy function first!")
  }
    
  # Checking modelname
  kaoliu_modelname <- tolower(kaoliu_modelname)
  kaoliu_modelname <- match.arg(kaoliu_modelname)
  model_modelname <- paste("model", kaoliu_modelname, sep = "_")
  
  dmunames <- datadea$dmunames
  nd <- length(dmunames) # number of dmus
  
  if (is.null(dmu_eval)) {
    dmu_eval <- 1:nd
  } else if (!all(dmu_eval %in% (1:nd))) {
    stop("Invalid set of DMUs to be evaluated (dmu_eval).")
  }
  names(dmu_eval) <- dmunames[dmu_eval]
  nde <- length(dmu_eval)
  
  arguments <- list(...)
  
  # Checking undesirable io and rts
  kkud1 <- FALSE # Flag for undesirable i/o basic non-directional
  kkud2 <- FALSE # Flag for undesirable i/o basic directional, RDM or SBM
  ud_inputs <- datadea$ud_inputs
  ud_outputs <- datadea$ud_outputs
  if (!is.null(ud_inputs) || !is.null(ud_outputs)) {
    if (kaoliu_modelname %in% c("basic", "fdh")) {
      if (is.null(arguments$orientation)) {
        arguments$orientation <- "io"
      }
      if (arguments$orientation != "dir") {
        kkud1 <- TRUE
        if (!is.null(ud_inputs) && (arguments$orientation != "io")) {
          warning("Undesirable (good) inputs with no input-oriented model.")
        }
        if (!is.null(ud_outputs) && (arguments$orientation != "oo")) {
          warning("Undesirable (bad) outputs with no output-oriented model.")
        }
        if (is.null(arguments$rts)) {
          arguments$rts <- "crs"
        }
        if (arguments$rts != "vrs") {
          #arguments$rts <- "vrs"
          warning("Returns to scale may be changed to variable (vrs) because there
                  is data with undesirable inputs/outputs.")
        }
      } else {
        kkud2 <- TRUE
      }
    } else if ((kaoliu_modelname == "sbmeff") || (kaoliu_modelname == "rdm")) {
      kkud2 <- TRUE
    } else {
      warning("This model does not take into account the undesirable feature for
              inputs/outputs.")
    }
  }
  ud_inputskkud2 <- NULL
  ud_outputskkud2 <- NULL
  if (kkud1) {
    datadea_old <- datadea
    res_und <- undesirable_basic(datadea = datadea, vtrans_i = arguments$vtrans_i,
                                 vtrans_o = arguments$vtrans_o)
    datadea <- res_und$u_datadea
    vtrans_i <- res_und$vtrans_i
    vtrans_o <- res_und$vtrans_o
  } else if (kkud2) {
    ud_inputskkud2 <- ud_inputs
    ud_outputskkud2 <- ud_outputs
  }
  
  input.mL <- datadea$input$mL
  output.mL <- datadea$output$mL
  input.mR <- datadea$input$mR
  output.mR <- datadea$output$mR
  input.dL <- datadea$input$dL
  output.dL <- datadea$output$dL
  input.dR <- datadea$input$dR
  output.dR <- datadea$output$dR
  inputnames <- rownames(input.mL)
  outputnames <- rownames(output.mL)
  ni <- length(inputnames) # number of inputs
  no <- length(outputnames) # number of outputs
  
  # Checking alpha
  if ((length(alpha) == 1) && (alpha > 1)) {
    alpha <- seq(from = 0, to = 1, length.out = alpha)
  } else if (any(alpha > 1) || any(alpha < 0)) {
    stop("Invalid alpha vector.")
  }
  nalpha <- length(alpha) # number of alpha-cuts
  alphacut <- vector(mode = "list", length = nalpha)
  names(alphacut) <- as.character(alpha)

  for (i in 1:nalpha) {

    # alpha-cuts
    
    a <- alpha[i]
    input.L <- input.mL - input.dL * (1 - a)
    input.U <- input.mR + input.dR * (1 - a)
    output.L <- output.mL - output.dL * (1 - a)
    output.U <- output.mR + output.dR * (1 - a)
    
    if (kkud2) { # Undesirable i/o for directional or SBM models
      if (!is.null(ud_inputs)){
        aux <- input.L[ud_inputs, ]
        input.L[ud_inputs, ] <- input.U[ud_inputs, ]
        input.U[ud_inputs, ] <- aux
      }
      if (!is.null(ud_outputs)){
        aux <- output.L[ud_outputs, ]
        output.L[ud_outputs, ] <- output.U[ud_outputs, ]
        output.U[ud_outputs, ] <- aux
      }
    }
    
    DMU.W <- vector(mode = "list", length = nde)
    names(DMU.W) <- dmunames[dmu_eval]
    DMU.B <- DMU.W
    
    for (j in 1:nde) {
      
      jj <- dmu_eval[j]
      
      # Worst data set for DMU jj
      input.W <- input.L
      input.W[, jj] <- input.U[, jj]
      output.W <- output.U
      output.W[, jj] <- output.L[, jj]
      datadea.W <- structure(list(input = input.W,
                                  output = output.W,
                                  dmunames = dmunames,
                                  nc_inputs = datadea$nc_inputs,
                                  nc_outputs = datadea$nc_outputs,
                                  nd_inputs = datadea$nd_inputs,
                                  nd_outputs = datadea$nd_outputs,
                                  ud_inputs = ud_inputskkud2,
                                  ud_outputs = ud_outputskkud2),
                             class = "deadata")
      deasol <- do.call(model_modelname, list(datadea = datadea.W, dmu_eval = jj, ...))
      DMU.W[[j]] <- deasol$DMU[[1]]
      
      # Best data set for DMU jj
      input.B <- input.U
      input.B[, jj] <- input.L[, jj]
      output.B <- output.L
      output.B[, jj] <- output.U[, jj]
      datadea.B <- structure(list(input = input.B,
                                  output = output.B,
                                  dmunames = dmunames,
                                  nc_inputs = datadea$nc_inputs,
                                  nc_outputs = datadea$nc_outputs,
                                  nd_inputs = datadea$nd_inputs,
                                  nd_outputs = datadea$nd_outputs,
                                  ud_inputs = ud_inputskkud2,
                                  ud_outputs = ud_outputskkud2),
                             class = "deadata")
      deasol <- do.call(model_modelname, list(datadea = datadea.B, dmu_eval = jj, ...))
      DMU.B[[j]] <- deasol$DMU[[1]]
      
    }
    
    if (kkud1) { 
      input.L <- datadea_old$input$mL - datadea_old$input$dL * (1 - a)
      input.U <- datadea_old$input$mR + datadea_old$input$dR * (1 - a)
      output.L <- datadea_old$output$mL - datadea_old$output$dL * (1 - a)
      output.U <- datadea_old$output$mR + datadea_old$output$dR * (1 - a)
      if (!is.null(DMU.W[[1]]$target_input)) {
        for (j in 1:nde) {
          DMU.W[[j]]$target_input[ud_inputs] <- vtrans_i - DMU.W[[j]]$target_input[ud_inputs]
          DMU.W[[j]]$target_output[ud_outputs] <- vtrans_o - DMU.W[[j]]$target_output[ud_outputs]
          DMU.B[[j]]$target_input[ud_inputs] <- vtrans_i - DMU.B[[j]]$target_input[ud_inputs]
          DMU.B[[j]]$target_output[ud_outputs] <- vtrans_o - DMU.B[[j]]$target_output[ud_outputs]
        }
      }
    }
    
    alphacut[[i]] <- list(input = list(Lower = input.L, Upper = input.U),
                          output = list(Lower = output.L, Upper = output.U),
                          DMU = list(Worst = DMU.W, Best = DMU.B))
    
  }
  
  deaOutput <- deasol
  
  if (kkud1) {
    datadea <- datadea_old
  }
  
  deaOutput$modelname <- paste("fuzzy_kaoliu", deasol$modelname, sep = "_")
  deaOutput$dmu_eval <- dmu_eval
  deaOutput$data <- datadea
  deaOutput$alpha <- alpha
  deaOutput$alphacut <- alphacut
  deaOutput$DMU <- NULL
  
  return(structure(deaOutput, class = "dea_fuzzy"))
  
}