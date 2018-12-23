#' @title Summary conventional DEA models.
#'   
#'   
#' @description Summary of the results obtained by a conventiona DEA model.
#' 
#' @param object An object of class \code{"dea"} obtained by a dea model function.
#' @param exportExcel Logical value. If TRUE (default) the results are also exported to an Excel file
#' @param filename Character string. Absolute filename (including path) of the exported Excel file. 
#'  If NULL, then the name of the file will be "ResultsDEA"+timestamp.xlsx.
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
#'                  dmus = 1, 
#'                  inputs = 2:6, 
#'                  outputs = 7:9 )
#' eval_pft <- model_basic(PFT, 
#'                         orientation = "io", 
#'                         rts = "crs")
#' summary(eval_pft, exporExcel = FALSE)
#' @references 
#' Charnes, A.; Cooper, W.W.; Rhodes, E. (1981). "Evaluating Program and Managerial Efficiency: An Application of Data Envelopment Analysis to Program Follow Through", Management Science, 27(6), 668-697. \url{https://pubsonline.informs.org/doi/abs/10.1287/mnsc.27.6.668}
#' @method summary dea
#' @import writexl
#' @importFrom dplyr summarise_at vars funs
#' @export
#' 

summary.dea <- function(object, exportExcel = TRUE, filename = NULL,...){
  
  if(!is.dea(object)){
    stop("Input should be of class dea!")
  }
 modelnames <-  c("basic", "additive", "addsupereff", "deaps", "fdh",
                  "multiplier", "nonradial","sbmeff", "sbmsupereff", 
                  "supereff","malmquist","cross_efficiency")
 modelname <- object$modelname
 # For CRAN - check pass
 Period <- vars <- ec <- mi <- mi <- funs <- DMU <- . <-  NULL
 if(!modelname %in% c("malmquist","cross_efficiency","bootstrap")){
 
 # Efficiencies
# if(!modelname %in% c("addsupereff")){
   eff <- efficiencies(object)
   eff <- data.frame(eff, stringsAsFactors = FALSE)
   eff <- data.frame(cbind(data.frame(DMU = rownames(eff)),eff), row.names = NULL)
# }else {
#   eff <- NULL
# }
 
 # slacks
 if(!modelname %in% c("multiplier")){
 s <- slacks(object) 
 s[sapply(s,is.null)] <- NULL
 s <- data.frame(s, stringsAsFactors = FALSE)
 s <- data.frame(cbind(data.frame(DMU = rownames(s)),s), row.names = NULL, stringsAsFactors = FALSE)
 
 } else{
   s <- NULL
 }
 # Lambdas
 lmbd <- lambdas(object)
 lamb <- data.frame(lmbd, stringsAsFactors = FALSE)
 lamb <- data.frame(cbind(data.frame(DMU = rownames(lamb)), lamb), 
                    row.names = NULL, stringsAsFactors = FALSE)
 
 # Targets
 tar <- targets(object)
 tar <- do.call(cbind,tar)
 #dimnames(tar)[[2]] <- paste("target",dimnames(tar)[[2]],sep = ".")
 tar <- data.frame(tar, stringsAsFactors = FALSE)
 tar <- data.frame(cbind(data.frame(DMU = rownames(tar)), tar), 
                   row.names = NULL, stringsAsFactors = FALSE)
 
 if(modelname == "multiplier"){
   mult <- multipliers(object)[1:2]
   mult <- do.call(cbind,mult)
  # dimnames(mult)[[2]] <- paste("multiplier",dimnames(mult)[[2]], sep = ".")
   mult <- data.frame(mult, stringsAsFactors = FALSE)
   mult <- data.frame(cbind(data.frame(DMU = object$data$dmunames), mult), 
                     row.names = NULL, stringsAsFactors = FALSE)
 }else{
   mult <- NULL
 }
 
 # References
 ref <- references(object) 
 #dmunames <- object$data$dmunames
 
 
 refnames <- unique(unlist(lapply(ref, function (x) names(x))))
 dmunames <- as.character(lamb$DMU)
 urefnames <- names(ref)
 
 
 RefMat <- matrix(0, nrow = length(dmunames), ncol = length(refnames),dimnames = list(dmunames,sort(refnames)))
 RefMat[urefnames,refnames] <- round(lmbd[urefnames, refnames],4)
 if(!modelname %in% c("addsupereff","sbmsupereff")){
   for(i in seq_along(refnames)){
     RefMat[refnames[i],refnames[i]] <- 1
   }
 }
 #refmat <- RefMat[urefnames,sort(refnames)]

 RefMatdf <- data.frame(cbind(data.frame(DMU = dmunames), data.frame(RefMat)), 
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
 returns <- data.frame(cbind(data.frame(DMU = rownames(returns)),returns), 
                       row.names = NULL, stringsAsFactors = FALSE)
 
 
 # Global data.frame
 
 
 dflist <- list(efficiencies = eff, 
                slacks = s,
                lambdas = lamb,
                targets = tar,
                multipliers = mult,
                returns = returns,
                references = RefMatdf) 
 dflist[sapply(dflist, is.null)] <- NULL
 
 if(exportExcel){
   if(is.null(filename)){
     filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
     filename <- gsub(" ","_",filename)
     filename <- gsub(":",".",filename)
   }
   write_xlsx(dflist, path = filename)
 }
 
 dflist <- lapply(dflist, function(x) x[-1])
 dffinal <- do.call(cbind,dflist)
 dffinal <- cbind(DMU = object$data$dmunames,dffinal)
 return(dffinal)
 }else if(modelname == "malmquist"){
   # Extract information about the data
   dmunames <- as.character(object$datadealist[[1]]$dmunames)
   periods <- names(object$datadealist)
   nper <- length(periods)
   
   
   # Create a list of data frames (each element is a period)
   df <- list()
   for (i in (1:(nper-1))){
     df[[i]]  <- data.frame(Period = periods[i+1], DMU = dmunames, ec = object$ec[i,], tc = object$tc[i,],pech = object$pech[i, ] ,
                            sech = object$sech[i,], mi = object$mi[i,])
   }
   # collapse the list into a data.frame
   dff <- do.call(rbind, df)
   rownames(dff) <- NULL

   # Geometric means by Period
   dff %>% group_by(Period) %>% summarise_at(vars(ec:mi),funs(geomean = exp(mean(log(.))))) %>% as.data.frame() -> dfsumPer
   colnames(dfsumPer) <- colnames(dff)[-2]
   # Geometric means by DMU
   dff %>% group_by(DMU) %>% summarise_at(vars(ec:mi),funs(geomean = exp(mean(log(.))))) %>% as.data.frame() -> dfsumDMU
   colnames(dfsumDMU) <- colnames(dff)[-1]
   
   res = list(Results = dff, means_by_period = dfsumPer, means_by_dmu = dfsumDMU)
   if(exportExcel){
     if(is.null(filename)){
       filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
       filename <- gsub(" ","_",filename)
       filename <- gsub(":",".",filename)
     }
     write_xlsx(res, path = filename)
   }
   return(res)
 }else if(modelname == "cross_efficiency"){
   nm <-  lapply(object, names)
   lst <- lapply(nm, function(x) "cross_eff" %in% x)
   lstce <- lst[sapply(lst,function(x) x)]
   dflist <- lapply(object[names(lstce)],function(x) x$cross_eff)
   
   dflist <- lapply(dflist, function(x) cbind(data.frame(DMU = dimnames(x)[[1]]),data.frame(x, row.names = NULL)))
   if(exportExcel){
     if(is.null(filename)){
       filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
       filename <- gsub(" ","_",filename)
       filename <- gsub(":",".",filename)
     }
     write_xlsx(dflist, path = filename)
   }
   return(dflist)
 } else{
   resMat <- cbind(object$score,object$score_bc,object$score-object$score_bc, object$CI)
   dimnames(resMat)[[2]] <- c("Score", "Bias-Corrected Score","Bias","CI Lower","CI Upper")
   resdf <- data.frame(cbind(data.frame(DMU = object$data$dmunames), data.frame(resMat)), 
                          row.names = NULL)
   if(exportExcel){
     if(is.null(filename)){
       filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
       filename <- gsub(" ","_",filename)
       filename <- gsub(":",".",filename)
     }
     write_xlsx(resdf, path = filename)
   }
   return(resdf)
 }
}