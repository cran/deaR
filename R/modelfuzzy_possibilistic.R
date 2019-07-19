#' @title Possibilistic Fuzzy DEA model.
#'   
#' @description Solve the possibilistic fuzzy DEA model proposed by León et. al (2003).
#' 
#' @usage modelfuzzy_possibilistic(datadea,
#'                          dmu_eval = NULL,
#'                          poss_modelname = c("basic"),
#'                          h = 1,
#'                          ...)
#' 
#' @param datadea The data, including DMUs, inputs and outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' @param poss_modelname a string containing the name of the model.
#' @param h A numeric vector with the h-levels (in [0,1]).
#' @param ... \code{dmu_ref}, \code{orientation}, \code{rts} and other model parameters.
#'   
#' @return An object of class \code{deadata_fuzzy}.
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
#' # Replication of results in Leon et. al (2003, p. 416)
#' data("Leon2003")
#' data_example <- read_data_fuzzy(Leon2003,
#'                                 inputs.mL = 2, 
#'                                 inputs.dL = 3, 
#'                                 outputs.mL = 4, 
#'                                 outputs.dL = 5)
#' result <- modelfuzzy_possibilistic(data_example, 
#'                                    h = seq(0, 1, by = 0.1), 
#'                                    orientation = "io", 
#'                                    rts = "vrs")
#' efficiencies(result)
#'  
#' @references
#' Emrouznejad, A.; Tavana, M.; Hatami-Marbini, A. (2014). “The State of the Art in Fuzzy Data Envelopment Analysis”, in A. Emrouznejad and M. Tavana (eds.), Performance Measurement with Fuzzy Data Envelopment Analysis. Studies in Fuzziness and Soft Computing 309. Springer, Berlin. \url{https://doi.org/10.1007/978-3-642-41372-8_1}
#' 
#' Hatami-Marbini, A.; Emrouznejad, A.; Tavana, M. (2011). "A Taxonomy and Review of the Fuzzy Data Envelopment Analysis Literature: Two Decades in the Making", European Journal of Operational Research, 214, 457–472. \url{https://doi.org/10.1016/j.ejor.2011.02.001}
#' 
#' Léon, T.; Liern, V. Ruiz, J.; Sirvent, I. (2003). "A Possibilistic Programming Approach to the Assessment of Efficiency with DEA Models", Fuzzy Sets and Systems, 139, 407–419. \url{https://doi.org/10.1016/S0165-0114(02)00608-5}
#' 
#' @seealso \code{\link{model_basic}}, \code{\link{modelfuzzy_kaoliu}}, \code{\link{modelfuzzy_guotanaka}}

#' 
#' @import lpSolve
#' 
#' @export
  
modelfuzzy_possibilistic <-
  function(datadea,
           dmu_eval = NULL,
           poss_modelname = c("basic"),
           h = 1,
           ...) {
 
  # Cheking whether datadea is of class "deadata_fuzzy" or not...  
  if (!is.deadata_fuzzy(datadea)) {
    stop("Data should be of class deadata_fuzzy. Run read_data_fuzzy function first!")
  }
    
  if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
    warning("This model does not take into account the undesirable feature for inputs/outputs.")
  }
  if (!is.null(datadea$nc_inputs) || !is.null(datadea$nc_outputs)) {
    warning("This model does not take into account the non-controllable feature for inputs/outputs.")
  }
  if (!is.null(datadea$nd_inputs) || !is.null(datadea$nd_outputs)) {
    warning("This model does not take into account the non-discretionary feature for inputs/outputs.")
  }
    
  # Checking modelname
  poss_modelname <- tolower(poss_modelname)
  poss_modelname <- match.arg(poss_modelname)
  model_modelname <- paste("model", poss_modelname, sep = "_")
  
  dmunames <- datadea$dmunames
  nd <- length(dmunames) # number of dmus
  
  if (is.null(dmu_eval)) {
    dmu_eval <- 1:nd
  } else if (all(dmu_eval %in% (1:nd)) == FALSE) {
    stop("Invalid set of DMUs to be evaluated (dmu_eval).")
  }
  names(dmu_eval) <- dmunames[dmu_eval]
  nde <- length(dmu_eval)
  
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
  
  datadea.mL <- structure(list(input = input.mL, output = output.mL, dmunames = dmunames), class = "deadata")
  datadea.mR <- structure(list(input = input.mR, output = output.mR, dmunames = dmunames), class = "deadata")
  
  # Checking h
  if (any(h > 1) || any(h < 0)){
    stop("Invalid h vector.")
  }
  h <- sort(unique(h))
  nh <- length(h) # number of h-levels
  hlevel <- vector(mode = "list", length = nh)
  names(hlevel) <- as.character(h)

  for (i in 1:nh) {

    # h-level
    
    a <- h[i]
    input.L <- input.mL - input.dL * (1 - a)
    input.U <- input.mR + input.dR * (1 - a)
    output.L <- output.mL - output.dL * (1 - a)
    output.U <- output.mR + output.dR * (1 - a)
    
    datadea.L <- structure(list(input = input.L, output = output.L, dmunames = dmunames), class = "deadata")
    datadea.U <- structure(list(input = input.U, output = output.U, dmunames = dmunames), class = "deadata")
    
    DMU <- vector(mode = "list", length = nde)
    names(DMU) <- dmunames[dmu_eval]
    
    for (j in 1:nde) {
      
      jj <- dmu_eval[j]
      
      lp.mL <- do.call(model_modelname, list(datadea = datadea.mL, dmu_eval = jj, returnlp = TRUE, ...))
      lp.mR <- do.call(model_modelname, list(datadea = datadea.mR, dmu_eval = jj, returnlp = TRUE, ...))
      lp.L <- do.call(model_modelname, list(datadea = datadea.L, dmu_eval = jj, returnlp = TRUE, ...))
      lp.U <- do.call(model_modelname, list(datadea = datadea.U, dmu_eval = jj, returnlp = TRUE, ...))
      
      obj <- lp.mL$DMU[[1]]$direction
      f.obj <- lp.mL$DMU[[1]]$objective.in
      
      f.con <- rbind(lp.mL$DMU[[1]]$const.mat,
                     lp.mR$DMU[[1]]$const.mat,
                     lp.L$DMU[[1]]$const.mat,
                     lp.U$DMU[[1]]$const.mat)
      
      f.dir <- c(lp.mL$DMU[[1]]$const.dir,
                 lp.mR$DMU[[1]]$const.dir,
                 lp.L$DMU[[1]]$const.dir,
                 lp.U$DMU[[1]]$const.dir)
      
      f.rhs <- c(lp.mL$DMU[[1]]$const.rhs,
                 lp.mR$DMU[[1]]$const.rhs,
                 lp.L$DMU[[1]]$const.rhs,
                 lp.U$DMU[[1]]$const.rhs)
      
      #if (poss_modelname == "fdh") {
      #  binary.vec = lp.mL$DMU[[1]]$binary.vec
      #} else {
      #  binary.vec = NULL
      #}
      #res <- lp(obj, f.obj, f.con, f.dir, f.rhs, binary.vec = binary.vec)
      
      res <- lp(obj, f.obj, f.con, f.dir, f.rhs)
      objval <- res$objval
      names(objval) <- "objval"
      res <- res$solution
      
      DMU[[j]] <- lp.mL$DMU[[1]]$var
      nv <- length(DMU[[j]])
      
      icontador <- 1
      for (iv in 1:nv) {
        nv2 <- length(DMU[[j]][[iv]])
        for (iv2 in 1:nv2) {
          DMU[[j]][[iv]][iv2] <- res[icontador]
          icontador <- icontador + 1
        }
      }
      
      #if (poss_modelname == "basic") {
      #}
      #if (poss_modelname == "additive") {
      #  DMU[[j]] <- c(objval, DMU[[j]])
      #} else if (poss_modelname == "nonradial"){
      #  mean_efficiency <- mean(DMU[[j]]$efficiency)
      #  DMU[[j]] <- c(mean_efficiency = mean_efficiency, DMU[[j]])
      #} else if (poss_modelname == "deaps") {
      #  mean_efficiency <- sum(lp.mL$weight * DMU[[j]]$efficiency) / sum(lp.mL$weight)
      #  DMU[[j]] <- c(mean_efficiency = mean_efficiency, DMU[[j]])
      #}
      
    }
    
    hlevel[[i]] <- list(input = list(Lower = input.L, Upper = input.U),
                        output = list(Lower = output.L, Upper = output.U),
                        DMU = DMU)
    
  }
    deaOutput <- list(modelname = paste("fuzzy_possibilistic", poss_modelname, sep = "_"),
                      orientation = lp.mL$orientation,
                      rts = lp.mL$rts,
                      L = lp.mL$L,
                      U = lp.mL$U,
                      h = h,
                      hlevel = hlevel,
                      data = datadea,
                      dmu_eval = dmu_eval,
                      dmu_ref = lp.mL$dmu_ref)
  #}
  
  return(structure(deaOutput, class = "dea_fuzzy"))
  
}