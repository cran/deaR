#' @title References
#'   
#' @description Extract the reference set for each DMU (inefficient DMUs and
#' efficicent DMUs that are combination of other efficient DMUs) from a DEA model solution.
#' 
#' @usage references(deasol,
#'            thr = 1e-4)
#' 
#' @param deasol Object of class \code{dea} obtained with some of the DEA model functions.
#' @param thr Tolerance threshold (for avoiding missdetection of  efficient 
#' DMUs due to roundoff errors)
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
#' # Replication results model DEA1 in Tomkins and Green (1988).
#' data("Departments")
#' # Calculate Total income
#' Departments$Total_income <- Departments[, 5] + Departments[, 6] + Departments[, 7] 
#' data_DEA1 <- read_data(Departments,
#'                        inputs = 9,
#'                        outputs = c(2, 3, 4, 12))
#' result <- model_basic(data_DEA1, 
#'                       orientation = "io",
#'                       rts = "crs")
#' references(result) # Table 3 (p.157)
#' 
#' @export

references <- function(deasol,
                       thr = 1e-4) {
  
  referencelist <- NULL
  
  if (is.dea(deasol)) {
    
    dmunames <- deasol$data$dmunames
    nd <- length(dmunames)
    dmu_eval <- deasol$dmu_eval
    nde <- length(dmu_eval)
    dmu_ref <- deasol$dmu_ref
    ndr <- length(dmu_ref)
    
    lamb <- lambdas(deasol)
    
    if(nrow(lamb) != nde) {
      dmu_eval <- which(sapply(deasol$DMU, function(x) !is.null(x$lambda)))
      nde <- length(dmu_eval)
    }
    
    # Alternatives in comment
    
    if (nde > 0) {
      
      #inputref <- matrix(deasol$data$input[, dmu_ref], nrow = ni) 
      #outputref <- matrix(deasol$data$output[, dmu_ref], nrow = no)
      #target_input <- t(inputref %*% t(lamb))
      #target_output <- t(outputref %*% t(lamb))
      testref <- rowSums(lamb) - rowSums(diag(nd)[dmu_eval, dmu_ref] * lamb) # is a DMU Inefficient?
      ineff_rows <- which(testref > thr)
      #ref <- lamb - diag(nd)[dmu_eval, dmu_ref]
      #colnames(ref) <- dmunames[dmu_ref] # psa
      #ineff_rows <- apply(X = ref, MARGIN = 1, FUN = function(x) any(x > thr)) # is a DMU Inefficient?
      #ineff_rows <- which(ineff_rows) 
      ineff_names <- dmunames[dmu_eval[ineff_rows]]
      nineff <- length(ineff_rows) # Number of ineff DMUs + eff DMUs that are combination of other eff DMUs
      
      referencelist <- vector(mode = "list", length = nineff) # Preallocate vector
      if (nineff == 0) {
        stop("There are no inefficient DMUs!")
      }
      for (i in 1:nineff) {
        aux <- which(lamb[ineff_rows[i], ] > thr)
        referencelist[[i]] <- lamb[ineff_rows[i], aux]
        
        names(referencelist[[i]]) <- names(dmu_ref)[aux]
        
      }
      
      names(referencelist) <- ineff_names
      
    }
    
  } else if (is.dea_fuzzy(deasol)) {
    
    dmunames <- deasol$data$dmunames
    nd <- length(dmunames)
    dmu_eval <- deasol$dmu_eval
    nde <- length(dmu_eval)
    dmu_ref <- deasol$dmu_ref
    ndr <- length(dmu_ref)
    dmunames_eval <- dmunames[dmu_eval]
    dmunames_ref <- dmunames[dmu_ref]
    ni <- length(deasol$data$input$mL[, 1])
    no <- length(deasol$data$output$mL[, 1])
    
    if (grepl("kaoliu", deasol$modelname)) {
      nalpha <- length(deasol$alpha)
      
      referencelist <- vector(mode = "list", length = nalpha) # Preallocate vector
      names(referencelist) <- names(deasol$alphacut)
      
      lamb <- lambdas(deasol)
      
      for (ii in 1:nalpha) {
        
        if(nrow(lamb$Worst[, , ii]) != nde) {
          dmu_eval <- which(sapply(deasol$alphacut[[ii]]$DMU$Worst, function(x) !is.null(x$lambda)))
          nde <- length(dmu_eval)
        }
        
        if (nde > 0) {
          
          ##### Worst
          
          testref <- rowSums(lamb$Worst[, , ii]) - rowSums(diag(nd)[dmu_eval, dmu_ref] * lamb$Worst[, , ii]) # is a DMU Inefficient?
          ineff_rows <- which(testref > thr)
          #ref <- lamb$Worst[, , ii] - diag(nd)[dmu_eval, dmu_ref]
          #colnames(ref) <- dmunames_ref # psa
          #ineff_rows <- apply(X = ref, MARGIN = 1, FUN = function(x) any(x > thr)) # is a DMU Inefficient?
          #ineff_rows <- which(ineff_rows) 
          ineff_names <- dmunames[dmu_eval[ineff_rows]]
          nineff <- length(ineff_rows) # Number of ineff DMUs + eff DMUs that are combination of other eff DMUs
          
          reflist <- vector(mode = "list", length = nineff) # Preallocate vector
          if (nineff == 0) {
            stop("There are no inefficient DMUs!")
          }
          for (i in 1:nineff) {
            aux <- which(lamb$Worst[ineff_rows[i], , ii] > thr)
            reflist[[i]] <- lamb$Worst[ineff_rows[i], aux, ii]
            
            names(reflist[[i]]) <- names(dmu_ref)[aux]
            
          }
          
          names(reflist) <- ineff_names
          
          referencelist[[ii]]$Worst <- reflist
          
          ##### Best
          
          if(nrow(lamb$Best[, , ii]) != nde) {
            dmu_eval <- which(sapply(deasol$alphacut[[ii]]$DMU$Best, function(x) !is.null(x$lambda)))
            nde <- length(dmu_eval)
          }
          
          if (nde > 0) {
            
            testref <- rowSums(lamb$Best[, , ii]) - rowSums(diag(nd)[dmu_eval, dmu_ref] * lamb$Best[, , ii]) # is a DMU Inefficient?
            ineff_rows <- which(testref > thr)
            #ref <- lamb$Best[, , ii] - diag(nd)[dmu_eval, dmu_ref]
            #colnames(ref) <- dmunames_ref # psa
            #ineff_rows <- apply(X = ref, MARGIN = 1, FUN = function(x) any(x > thr)) # is a DMU Inefficient?
            #ineff_rows <- which(ineff_rows) 
            ineff_names <- dmunames[dmu_eval[ineff_rows]]
            nineff <- length(ineff_rows) # Number of ineff DMUs + eff DMUs that are combination of other eff DMUs
            
            reflist <- vector(mode = "list", length = nineff) # Preallocate vector
            if (nineff == 0) {
              stop("There are no inefficient DMUs!")
            }
            for (i in 1:nineff) {
              aux <- which(lamb$Best[ineff_rows[i], , ii] > thr)
              reflist[[i]] <- lamb$Best[ineff_rows[i], aux, ii]
              
              names(reflist[[i]]) <- names(dmu_ref)[aux]
              
            }
            
            names(reflist) <- ineff_names
            
            referencelist[[ii]]$Best <- reflist
        
          }
          
        }
        
      }
      
    } else if (grepl("possibilistic", deasol$modelname)) {
      nh <- length(deasol$h)
      
      referencelist <- vector(mode = "list", length = nh) # Preallocate vector
      names(referencelist) <- names(deasol$hlevel)
      
      lamb <- lambdas(deasol)
      
      for (ii in 1:nh) {
        
        if(nrow(lamb[, , ii]) != nde) {
          dmu_eval <- which(sapply(deasol$hlevel[[ii]]$DMU, function(x) !is.null(x$lambda)))
          nde <- length(dmu_eval)
        }
        
        if (nde > 0) {
          
          testref <- rowSums(lamb[, , ii]) - rowSums(diag(nd)[dmu_eval, dmu_ref] * lamb[, , ii]) # is a DMU Inefficient?
          ineff_rows <- which(testref > thr)
          #ref <- lamb[, , ii] - diag(nd)[dmu_eval, dmu_ref]
          #colnames(ref) <- dmunames_ref # psa
          #ineff_rows <- apply(X = ref, MARGIN = 1, FUN = function(x) any(x > thr)) # is a DMU Inefficient?
          #ineff_rows <- which(ineff_rows) 
          ineff_names <- dmunames[dmu_eval[ineff_rows]]
          nineff <- length(ineff_rows) # Number of ineff DMUs + eff DMUs that are combination of other eff DMUs
          
          reflist <- vector(mode = "list", length = nineff) # Preallocate vector
          if (nineff == 0) {
            stop("There are no inefficient DMUs!")
          }
          for (i in 1:nineff) {
            aux <- which(lamb[ineff_rows[i], , ii] > thr)
            reflist[[i]] <- lamb[ineff_rows[i], aux, ii]
            
            names(reflist[[i]]) <- names(dmu_ref)[aux]
            
          }
          
          names(reflist) <- ineff_names
          
          referencelist[[ii]] <- reflist
          
        }
        
      }
      
    } else {
      stop("No target parameters in this solution!")
    }
    
    
  } else {
    
    stop("Input should be a dea or dea_fuzzy class object!")
    
  }
  
  return(referencelist)
}