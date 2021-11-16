#' @title Summary conventional DEA models.
#'   
#'   
#' @description Summary of the results obtained by a conventiona DEA model.
#' 
#' 
#' @param object An object of class \code{"dea"} obtained by a dea model function.
#' @param exportExcel Logical value. If TRUE (default) the results are also exported to an Excel file
#' @param filename Character string. Absolute filename (including path) of the exported Excel file. 
#'  If NULL, then the name of the file will be "ResultsDEA"+timestamp.xlsx.
#' @param returnList Logical value. If TRUE then the results are given as a list of data frames. 
#'  If FALSE (default) all the data frames are merged into a single data frame.
#' @param ... Ignored. Used for compatibility issues.
#' 
#'   
#' @return Depending on the model it returns a single data.frame containing: efficiencies, 
#' slacks, lambdas, targets, references or a list of data.frames with the cross-efficiencies computed 
#' with different methods (Arbitrary, Method II or Method III (see CITA)) or, in case the model is a
#'  malmquist index, a single data.frame with the coefficients for the different periods.       
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
#' data("PFT1981") 
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ] 
#' PFT <- read_data(PFT, 
#'                  inputs = 2:6, 
#'                  outputs = 7:9 )
#' eval_pft <- model_basic(PFT, 
#'                         orientation = "io", 
#'                         rts = "crs")
#' summary(eval_pft, exportExcel = FALSE)
#' @references 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1981). "Evaluating Program and Managerial 
#' Efficiency: An Application of Data Envelopment Analysis to Program Follow Through", 
#' Management Science, 27(6), 668-697. 
#' \url{https://pubsonline.informs.org/doi/abs/10.1287/mnsc.27.6.668}
#' @method summary dea
#' @import writexl
#' @importFrom dplyr summarise_at vars funs
#' @export
#' 

summary.dea <- function(object, 
                        exportExcel = TRUE,
                        filename = NULL, 
                        returnList = FALSE,
                        ...) {
  
  if (!is.dea(object)) {
    stop("Input should be of class dea!")
  }
  modelnames <-
    c(
      "basic",
      "additive",
      "addsupereff",
      "deaps",
      "fdh",
      "multiplier",
      "nonradial",
      "sbmeff",
      "sbmsupereff",
      "supereff",
      "malmquist",
      "cross_efficiency",
      "profit"
    )
  modelname <- object$modelname
  # For CRAN - check pass
  Period <- vars <- ec <- mi <- mi <- funs <- DMU <- . <-  NULL
  
  if (!modelname %in% c("malmquist", "cross_efficiency", "bootstrap", "profit")) {
    # All models except malmquist, ce, bootstrap and profit -------
    # Efficiencies
    # if(!modelname %in% c("addsupereff")){
    eff <- efficiencies(object)
    eff <- data.frame(eff, stringsAsFactors = FALSE)
    eff <-
      data.frame(cbind(data.frame(DMU = rownames(eff)), eff), row.names = NULL)
    # }else {
    #   eff <- NULL
    # }
    
    # slacks
    if (!modelname %in% c("multiplier")) {
      s <- slacks(object)
      s[sapply(s, is.null)] <- NULL
      s <- data.frame(s, stringsAsFactors = FALSE)
      s <-
        data.frame(cbind(data.frame(DMU = rownames(s)), s),
                   row.names = NULL,
                   stringsAsFactors = FALSE)
      
    } else {
      s <- NULL
    }
    # Lambdas
    lmbd <- lambdas(object)
    lamb <- data.frame(lmbd, stringsAsFactors = FALSE)
    lamb <- data.frame(cbind(data.frame(DMU = rownames(lamb)), lamb),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
    
    # Targets
    tar <- targets(object)
    tar <- do.call(cbind, tar)
    #dimnames(tar)[[2]] <- paste("target",dimnames(tar)[[2]],sep = ".")
    tar <- data.frame(tar, stringsAsFactors = FALSE)
    tar <- data.frame(cbind(data.frame(DMU = rownames(tar)), tar),
                      row.names = NULL,
                      stringsAsFactors = FALSE)
    
    if (modelname == "multiplier") {
      mult <- multipliers(object)[1:2]
      mult <- do.call(cbind, mult)
      # dimnames(mult)[[2]] <- paste("multiplier",dimnames(mult)[[2]], sep = ".")
      mult <- data.frame(mult, stringsAsFactors = FALSE)
      mult <-
        data.frame(cbind(data.frame(DMU = object$data$dmunames), mult),
                   row.names = NULL,
                   stringsAsFactors = FALSE)
    } else {
      mult <- NULL
    }
    
    # References
    ref <- references(object)
    #dmunames <- object$data$dmunames
    
    
    refnames <- unique(unlist(lapply(ref, function (x)
      names(x))))
    dmunames <- as.character(lamb$DMU)
    urefnames <- names(ref)
    
    
    RefMat <-
      matrix(
        0,
        nrow = length(dmunames),
        ncol = length(refnames),
        dimnames = list(dmunames, sort(refnames))
      )
    RefMat[urefnames, refnames] <- round(lmbd[urefnames, refnames], 4)
    if (!modelname %in% c("addsupereff", "sbmsupereff")) {
      for (i in seq_along(refnames)) {
        if (refnames[i] %in% dimnames(RefMat)[[1]] & refnames[i] %in% dimnames(RefMat)[[2]]){
          RefMat[refnames[i], refnames[i]] <- 1
        }
      }
    }
    #refmat <- RefMat[urefnames,sort(refnames)]
    
    RefMatdf <-
      data.frame(cbind(data.frame(DMU = dmunames), data.frame(RefMat)),
                 row.names = NULL)
    
    #
    # refmat  <- matrix(0, nrow = length(ref),
    #                   ncol = length(refnames),
    #                   dimnames = list(names(ref), refnames))
    # refmat[names(ref),refnames] <- round(lmbd[names(ref), refnames],4)
    # refmat <- refmat[,sort(dimnames(refmat)[[2]])]
    
    
    
    # Returns
    returns <- rts(object)
    returns <- data.frame(returns)
    returns <-
      data.frame(cbind(data.frame(DMU = rownames(returns)), returns),
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    
    
    # Global data.frame
    
    
    dflist <- list(
      efficiencies = eff,
      slacks = s,
      lambdas = lamb,
      targets = tar,
      multipliers = mult,
      returns = returns,
      references = RefMatdf
    )
    dflist[sapply(dflist, is.null)] <- NULL
    
    if (exportExcel) {
      if (is.null(filename)) {
        filename <- paste("ResultsDEA", Sys.time(), ".xlsx", sep = "")
        filename <- gsub(" ", "_", filename)
        filename <- gsub(":", ".", filename)
      }
      write_xlsx(dflist, path = filename)
    }
    
    #dflist <- lapply(dflist, function(x) x[-1])
    if (returnList) {
      return(dflist)
    } else {
      dffinal <- do.call(cbind, dflist)
      dffinal <- cbind(DMU = names(object$dmu_eval), dffinal)
      return(dffinal)
    }
  } else if (modelname == "malmquist") {
    # Malmquist model -----
    # Extract information about the data
    #dmunames <- as.character(object$datadealist[[1]]$dmunames)
    dmunames <- names(object$dmu_eval)
    periods <- names(object$datadealist)
    nper <- length(periods)
    
    
    # Create a list of data frames (each element is a period)
    df <- list()
    reslist <- object[1:8]
    reslist[sapply(reslist,is.null)] <- NULL
    
    for (i in (1:(nper - 1))) {
      df[[i]]  <-
        cbind(data.frame(
          Period = periods[i + 1],
          DMU = dmunames),
          sapply(reslist, function(x) x[i,]))
        # data.frame(
        #   Period = periods[i + 1],
        #   DMU = dmunames,
        #   ec = object$ec[i, ],
        #   tc = object$tc[i, ],
        #   pech = object$pech[i,] ,
        #   sech = object$sech[i, ],
        #   mi = object$mi[i, ]
        # )
    }
    # collapse the list into a data.frame
    dff <- do.call(rbind, df)
    rownames(dff) <- NULL
    cnames <- colnames(dff)
    # Geometric means by Period vars(3:ncol(dff))
    dff %>% group_by(Period) %>% summarise_at(vars(cnames[3]:cnames[ncol(dff)]),
                                              list(geomean = ~exp(mean(log(
      .
    ))))) %>% as.data.frame() -> dfsumPer
    colnames(dfsumPer) <- colnames(dff)[-2]
    # Geometric means by DMU
    dff %>% group_by(DMU) %>% summarise_at(vars(cnames[3]:cnames[ncol(dff)]), 
                                           list(geomean = ~exp(mean(log(
      .
    ))))) %>% as.data.frame() -> dfsumDMU
    colnames(dfsumDMU) <- colnames(dff)[-1]
    
    res = list(
      Results = dff,
      means_by_period = dfsumPer,
      means_by_dmu = dfsumDMU
    )
    if (exportExcel) {
      if (is.null(filename)) {
        filename <- paste("ResultsDEA", Sys.time(), ".xlsx", sep = "")
        filename <- gsub(" ", "_", filename)
        filename <- gsub(":", ".", filename)
      }
      write_xlsx(res, path = filename)
    }
    return(res)
  } else if (modelname == "cross_efficiency") {
    # Cross - efficiency -----
    nm <-  lapply(object, names)
    lst <- lapply(nm, function(x)
      "cross_eff" %in% x)
    lstce <- lst[sapply(lst, function(x)
      x)]
    dflist <- lapply(object[names(lstce)], function(x)
      x$cross_eff)
    
    dflist <-
      lapply(dflist, function(x)
        cbind(data.frame(DMU = dimnames(x)[[1]]), data.frame(x, row.names = NULL)))
    if (exportExcel) {
      if (is.null(filename)) {
        filename <- paste("ResultsDEA", Sys.time(), ".xlsx", sep = "")
        filename <- gsub(" ", "_", filename)
        filename <- gsub(":", ".", filename)
      }
      write_xlsx(dflist, path = filename)
    }
    return(dflist)
  } else if (modelname == "profit") {
    # Profit model -------
    modeltype <-
      ifelse(
        !is.null(object$price_input),
        ifelse(
          is.null(object$price_output),
          "price_input",
          "price_input_output"
        ),
        "price_output"
      )
    
    switch (
      modeltype,
      price_input = {
        effname <- "cost_efficiency"
        objname <- "minimum_cost"
      },
      price_output = {
        effname <- "revenue_efficiency"
        objname <- "maximum_revenue"
      },
      price_input_output = {
        effname <- "profit_efficiency"
        objname <- "maximum_profit"
      }
    )
    # Efficiencies
    eff <- efficiencies(object)
    eff <- data.frame(eff, stringsAsFactors = FALSE)
    colnames(eff) <- effname
    eff <-
      data.frame(cbind(data.frame(DMU = rownames(eff)), eff), row.names = NULL)
    # Lambdas
    lmbd <- lambdas(object)
    lamb <- data.frame(lmbd, stringsAsFactors = FALSE)
    lamb <-
      data.frame(cbind(data.frame(DMU = rownames(lamb)), lamb),
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    # Objective value
    objval <- unlist(lapply(object$DMU, function(x)
      x$objval))
    objval <- data.frame(objval, stringsAsFactors = FALSE)
    colnames(objval) <- objname
    objval <-
      data.frame(cbind(data.frame(DMU = rownames(objval)), objval), row.names = NULL)
    # RTS
    returns <- rts(object)
    returns <- data.frame(returns)
    returns <-
      data.frame(cbind(data.frame(DMU = rownames(returns)), returns),
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    # references
    ref <- references(object)
    refnames <- unique(unlist(lapply(ref, function (x)
      names(x))))
    dmunames <- as.character(lamb$DMU)
    urefnames <- names(ref)
    
    
    
    RefMat <-
      matrix(
        0,
        nrow = length(dmunames),
        ncol = length(refnames),
        dimnames = list(dmunames, sort(refnames))
      )
    RefMat[urefnames, refnames] <- round(lmbd[urefnames, refnames], 4)
    
    RefMatdf <-
      data.frame(cbind(data.frame(DMU = dmunames), data.frame(RefMat)),
                 row.names = NULL)
    
    
    # Optimal i/o
    switch(
      modeltype,
      price_input = {
        optimio <-
          do.call(rbind, lapply(object$DMU, function(x)
            x$optimal_input))
      },
      price_output = {
        optimio <-
          do.call(rbind, lapply(object$DMU, function(x)
            x$optimal_input))
      },
      price_input_output = {
        oi <- do.call(rbind, lapply(object$DMU, function(x)
          x$optimal_input))
        oo <-
          do.call(rbind, lapply(object$DMU, function(x)
            x$optimal_output))
        optimio <- cbind(oi, oo)
      }
    )
    optimio <-
      data.frame(cbind(data.frame(DMU = rownames(optimio)), optimio),
                 row.names = NULL,
                 stringsAsFactors = FALSE)
    
    
    dflist <- list(
      efficiencies = eff,
      objval = objval,
      lambdas = lamb,
      returns = returns,
      optimio = optimio,
      references = RefMatdf
    )
    dflist[sapply(dflist, is.null)] <- NULL
    if (exportExcel) {
      if (is.null(filename)) {
        filename <- paste("ResultsDEA", Sys.time(), ".xlsx", sep = "")
        filename <- gsub(" ", "_", filename)
        filename <- gsub(":", ".", filename)
      }
      write_xlsx(dflist, path = filename)
    }
    
   # dflist <- lapply(dflist, function(x)
    #  x[-1])
    if (returnList) {
      return(dflist)
    } else {
      dffinal <- do.call(cbind, dflist)
      dffinal <- cbind(DMU = object$data$dmunames, dffinal)
      return(dffinal)
    }
  } else {
    # Bootstrap -----
    resMat <-
      cbind(object$score,
            object$score_bc,
            object$score - object$score_bc,
            object$CI)
    dimnames(resMat)[[2]] <-
      c("Score", "Bias-Corrected Score", "Bias", "CI Lower", "CI Upper")
    resdf <-
      data.frame(cbind(data.frame(DMU = object$data$dmunames), data.frame(resMat)),
                 row.names = NULL)
    if (exportExcel) {
      if (is.null(filename)) {
        filename <- paste("ResultsDEA", Sys.time(), ".xlsx", sep = "")
        filename <- gsub(" ", "_", filename)
        filename <- gsub(":", ".", filename)
      }
      write_xlsx(resdf, path = filename)
    }
    return(resdf)
  }
}
