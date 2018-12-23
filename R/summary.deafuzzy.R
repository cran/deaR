#' @title Summary Fuzzy DEA models.
#'   
#'   
#' @description Summary of the results obtained by a fuzzy DEA model.
#' @param object An object of class \code{"dea_fuzzy"} obtained with a fuzzy dea 
#'  model function (\code{modelfuzzy_guotanaka}, \code{modelfuzzy_kaoliu},
#'  \code{modelfuzzy_possibilistic}).
#' @param ... Extra options
#' @param exportExcel Logical value. If TRUE (default) the results are also exported to an Excel file
#' @param filename Character string. Absolute filename (including path) of the exported Excel file. 
#'  If NULL, then the name of the file will be "ResultsDEA"+timestamp.xlsx.
#' 
#'   
#' @return If the model is that from Guo and Tanaka (\code{modelfuzzy_guotanaka}), it returns a data.frame
#' with columns: DMU, alpha cuts and efficiencies.
#' For the possibilistic model (\code{modelfuzzy_possibilistic}) it returns a data.frame with columns:
#' DMU, alpha-cuts, efficiencies and the corresponding lambda values
#' For the Kao and Liu model (\code{modelfuzzy_kaoliu}), the result may depend on the crisp sub-model used. 
#' It will contain a data.frame with the efficiencies (if any), the slacks and superslacks (if any), 
#' the lambda values and the targets.
#' 
#' If \code{exportExcel} is TRUE, then an Excel file will be created containing as many sheets as necessary 
#' depending on the variables returned.
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
#' @examples
#' data("Leon2003")
#' data_example <- read_data_fuzzy(Leon2003,
#'                                 dmus = 1, 
#'                                 inputs.mL = 2, 
#'                                 inputs.dL = 3, 
#'                                 outputs.mL = 4, 
#'                                 outputs.dL = 5)
#' result <- modelfuzzy_possibilistic(data_example, 
#'                                    h = seq(0, 1, by = 0.1), 
#'                                    orientation = "io", 
#'                                    rts = "vrs")
#' summary(result)
#' @references 
#' León, T.; Liern, V. Ruiz, J.; Sirvent, I. (2003). "A Possibilistic Programming Approach to the Assessment of Efficiency with DEA Models", Fuzzy Sets and Systems, 139, 407–419. \url{https://doi.org/10.1016/S0165-0114(02)00608-5}
#' @method summary dea_fuzzy
#' @import writexl tidyr
#' @export
#' 

summary.dea_fuzzy <- function(object, ..., exportExcel = TRUE, filename = NULL){
  
  if(!is.dea_fuzzy(object)){
    stop("Input should be of class dea_fuzzy!")
  }
 
 modelname <- object$modelname
 # For CRAN - check pass
 DMU <- NULL
 # Guo-Tanaka -----------------------
 if(modelname == "fuzzy_guotanaka"){
   eff <- efficiencies(object)
   effmat <- do.call(rbind,lapply(seq(dim(eff)[3]), function(x) eff[,,x] ))
   effdf <- cbind(data.frame(DMU = dimnames(effmat)[[1]],
                             hlevel = rep(object$h,each = dim(eff)[1])), 
                  data.frame(effmat, row.names = NULL))
   
   if(exportExcel){
     if(is.null(filename)){
       filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
       filename <- gsub(" ","_",filename)
       filename <- gsub(":",".",filename)
     }
     write_xlsx(effdf, path = filename)
   }
   return(effdf)
   # Possibilistic --------------------
 }else if(modelname == "fuzzy_possibilistic_basic"){
   # Efficiencies
   eff <- efficiencies(object)
   eff <- cbind(data.frame(DMU = dimnames(eff)[[1]]), 
                data.frame(eff, row.names = NULL))
   eff %>% gather(key = "hlevel", value = "efficiency",-DMU) -> eff
   eff$hlevel <- rep(object$h, each = length(object$data$dmunames))
   
   # Lambdas
   lamb <- lambdas(object)
   lamblist <- lapply(seq(dim(lamb)[3]), function(x) lamb[,,x])
   lambmat <- do.call(rbind,lamblist)
   
   df <- cbind(eff, data.frame(lambmat, row.names = NULL))
   if(exportExcel){
     if(is.null(filename)){
       filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
       filename <- gsub(" ","_",filename)
       filename <- gsub(":",".",filename)
     }
     write_xlsx(df, path = filename)
   }
   return(df)
   # Kao - Liu ---------------------
 }else{
   modelkl <- strsplit(object$modelname,"_")[[1]][3]
   
   # Efficiencies ----------
   if(!modelkl %in% c("additive","addsupereff")){
     eff <- efficiencies(object)
     
     # Radial models ----------------
     if(!modelkl %in% c("nonradial","deaps")){
       
       eff.Worst <- data.frame(eff$Worst, stringsAsFactors = FALSE)
       eff.Worst <- data.frame(cbind(data.frame(DMU = rownames(eff.Worst)),
                                     eff.Worst), 
                               row.names = NULL)
       eff.Worst %>% gather(key = "alphacut", 
                            value = "efficiency.Worst",-DMU) -> eff.Worst
       eff.Worst$alphacut <- rep(object$alpha, 
                                  each = length(object$data$dmunames))
       
       eff.Best <- data.frame(eff$Best, stringsAsFactors = FALSE)
       eff.Best <- data.frame(cbind(data.frame(DMU = rownames(eff.Best)),eff.Best), row.names = NULL)
       eff.Best %>% gather(key = "alphacut", value = "efficiency.Best",-DMU) -> eff.Best
       eff.Best$alphacut <- rep(object$alpha, each = length(object$data$dmunames))
       eff.df <- merge(eff.Worst,eff.Best, by = c("DMU","alphacut"))
     }else{
       # Non - radial models-----------------
       effmat.Worst <- do.call(rbind,lapply(seq(dim(eff$Worst)[3]), function(x) eff$Worst[,,x] ))
       effdf.Worst <- cbind(data.frame(DMU = dimnames(effmat.Worst)[[1]], 
                                   alphacut = rep(object$alpha,each = dim(eff$Worst)[1])), 
                        data.frame(effmat.Worst, row.names = NULL))
       colnames(effdf.Worst)[3:(ncol(effdf.Worst))] <- paste("eff",colnames(effdf.Worst)[3:(ncol(effdf.Worst))],"Worst",sep = ".")
       
       effmat.Best <- do.call(rbind,lapply(seq(dim(eff$Best)[3]), function(x) eff$Best[,,x] ))
       effdf.Best <- cbind(data.frame(DMU = dimnames(effmat.Best)[[1]], 
                                   alphacut = rep(object$alpha,each = dim(eff$Best)[1])), 
                        data.frame(effmat.Best, row.names = NULL))
       colnames(effdf.Best)[3:(ncol(effdf.Best))] <- paste("eff",colnames(effdf.Best)[3:(ncol(effdf.Best))],"Best",sep = ".")
       eff.df  <- merge(effdf.Worst,effdf.Best, by = c("alphacut","DMU"))
       # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
       srtidx <- (3:ncol(eff.df))
       srtidx <- t(matrix(srtidx, ncol = 2))
       dim(srtidx) <- c(1,length(srtidx))
       eff.df <- eff.df[,c(2,1,srtidx)]
       
     }
   }else {
     eff.df <- NULL
   }
   
   # Slacks -------------------------
     s <- slacks(object)
     s[sapply(s, is.null)] <- NULL # Remove NULL fields from slacks list
     dmunames <- object$data$dmunames
     # Case 1: both input and output slacks are present -----------------
     
     if(!modelkl %in% c("nonradial","deaps")){
    
     s.i.Worst <- do.call(rbind,lapply(seq(dim(s$slack_input.W)[3]), 
                                       function(x) matrix(s$slack_input.W[,,x],
                                                          nrow = length(dmunames),
                                                          dimnames = list(dmunames,dimnames(s$slack_input.W)[[2]]))))
     s.i.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.i.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(s$slack_input.W)[1])), 
                          data.frame(s.i.Worst, row.names = NULL))
     colnames(s.i.Worst)[3:(ncol(s.i.Worst))] <- paste("slack",
                                                       colnames(s.i.Worst)[3:(ncol(s.i.Worst))],
                                                       "Worst",sep = ".")
     
     s.o.Worst <- do.call(rbind,lapply(seq(dim(s$slack_output.W)[3]), 
                                       function(x) matrix(s$slack_output.W[,,x],
                                                          nrow = length(dmunames),
                                                          dimnames = list(dmunames,dimnames(s$slack_output.W)[[2]]))))
     s.o.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.o.Worst)[[1]], 
                                   alphacut = rep(object$alpha,each = dim(s$slack_output.W)[1])), 
                        data.frame(s.o.Worst, row.names = NULL))
     colnames(s.o.Worst)[3:(ncol(s.o.Worst))] <- paste("slack",
                                                       colnames(s.o.Worst)[3:(ncol(s.o.Worst))],
                                                       "Worst",sep = ".")
     
     s.i.Best <- do.call(rbind,lapply(seq(dim(s$slack_input.B)[3]), 
                                      function(x) matrix(s$slack_input.B[,,x],
                                                         nrow = length(dmunames),
                                                         dimnames = list(dmunames,dimnames(s$slack_input.B)[[2]])) ))
     s.i.Best <- cbind(data.frame(DMU = dimnames(s.i.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(s$slack_input.B)[1])), 
                         data.frame(s.i.Best, row.names = NULL))
     colnames(s.i.Best)[3:(ncol(s.i.Best))] <- paste("slack",
                                                     colnames(s.i.Best)[3:(ncol(s.i.Best))],
                                                     "Best",sep = ".")
     s.o.Best <- do.call(rbind,lapply(seq(dim(s$slack_output.B)[3]), 
                                      function(x) matrix(s$slack_output.B[,,x],
                                                         nrow = length(dmunames),
                                                         dimnames = list(dmunames,dimnames(s$slack_output.B)[[2]])) ))
     s.o.Best <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.o.Best)[[1]], 
                                  alphacut = rep(object$alpha,each = dim(s$slack_output.B)[1])), 
                       data.frame(s.o.Best, row.names = NULL))
     colnames(s.o.Best)[3:(ncol(s.o.Best))] <- paste("slack",
                                                     colnames(s.o.Best)[3:(ncol(s.o.Best))],
                                                     "Best",sep = ".")
     
     s.i.df  <- merge(s.i.Worst,s.i.Best, by = c("alphacut","DMU"))
     # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
     srtidx <- (3:ncol(s.i.df))
     srtidx <- t(matrix(srtidx, ncol = 2))
     dim(srtidx) <- c(1,length(srtidx))
     s.i.df <- s.i.df[,c(2,1,srtidx)]
     
     s.o.df  <- merge(s.o.Worst,s.o.Best, by = c("alphacut","DMU"))
     # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
     srtidx <- (3:ncol(s.o.df))
     srtidx <- t(matrix(srtidx, ncol = 2))
     dim(srtidx) <- c(1,length(srtidx))
     s.o.df <- s.o.df[,c(2,1,srtidx)]
     s.df <- cbind(s.i.df, s.o.df[,3:ncol(s.o.df)])
     
     
     
   }else{
     # Case 2: only either input or output slacks are present (but not both) --------------
     if(object$orientation == "io"){
       s.o.Worst <- do.call(rbind,lapply(seq(dim(s$slack_output.W)[3]), 
                                         function(x) matrix(s$slack_output.W[,,x],
                                                            nrow = length(dmunames),
                                                            dimnames = list(dmunames,dimnames(s$slack_output.W)[[2]]))))
       s.o.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.o.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(s$slack_output.W)[1])), 
                          data.frame(s.o.Worst, row.names = NULL))
       colnames(s.o.Worst)[3:(ncol(s.o.Worst))] <- paste("slack",
                                                         colnames(s.o.Worst)[3:(ncol(s.o.Worst))],
                                                         "Worst",sep = ".")
       
       
       
       s.o.Best <- do.call(rbind,lapply(seq(dim(s$slack_output.B)[3]), 
                                        function(x) matrix(s$slack_output.B[,,x],
                                                           nrow = length(dmunames),
                                                           dimnames = list(dmunames,dimnames(s$slack_output.B)[[2]])) ))
       s.o.Best <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.o.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(s$slack_output.B)[1])), 
                         data.frame(s.o.Best, row.names = NULL))
       colnames(s.o.Best)[3:(ncol(s.o.Best))] <- paste("slack",
                                                       colnames(s.o.Best)[3:(ncol(s.o.Best))],
                                                       "Best",sep = ".")
       
       s.o.df  <- merge(s.o.Worst,s.o.Best, by = c("alphacut","DMU"))
       # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
       srtidx <- (3:ncol(s.o.df))
       srtidx <- t(matrix(srtidx, ncol = 2))
       dim(srtidx) <- c(1,length(srtidx))
       s.df <- s.o.df[,c(2,1,srtidx)]
     }else {
       s.i.Worst <- do.call(rbind,lapply(seq(dim(s$slack_input.Wslack_input.W)[3]), 
                                         function(x) matrix(s$slack_input.W[,,x],
                                                            nrow = length(dmunames),
                                                            dimnames = list(dmunames,dimnames(s$slack_input.W)[[2]]))))
       s.i.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(s.i.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(s$slack_input.W)[1])), 
                          data.frame(s.i.Worst, row.names = NULL))
       colnames(s.i.Worst)[3:(ncol(s.i.Worst))] <- paste("slack",
                                                         colnames(s.i.Worst)[3:(ncol(s.i.Worst))],
                                                         "Worst",sep = ".")
       
       s.i.Best <- do.call(rbind,lapply(seq(dim(s$slack_input.B)[3]), 
                                        function(x) matrix(s$slack_input.B[,,x],
                                                           nrow = length(dmunames),
                                                           dimnames = list(dmunames,dimnames(s$slack_input.B)[[2]])) ))
       s.i.Best <- cbind(data.frame(DMU = dimnames(s.i.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(s$slack_input.B)[1])), 
                         data.frame(s.i.Best, row.names = NULL))
       colnames(s.i.Best)[3:(ncol(s.i.Best))] <- paste("slack",
                                                       colnames(s.i.Best)[3:(ncol(s.i.Best))],
                                                       "Best",sep = ".")
       s.i.df  <- merge(s.i.Worst,s.i.Best, by = c("alphacut","DMU"))
       # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
       srtidx <- (3:ncol(s.i.df))
       srtidx <- t(matrix(srtidx, ncol = 2))
       dim(srtidx) <- c(1,length(srtidx))
       s.df <- s.i.df[,c(2,1,srtidx)]
     }
   }
     
     # Case 3: extra option for supers! Superslacks present -------------
     if(modelkl %in% c("addsupereff","sbmsupereff")){
       supers.i.Worst <- do.call(rbind,lapply(seq(dim(s$superslack_input.W)[3]), 
                                         function(x) matrix(s$superslack_input.W[,,x],
                                                            nrow = length(dmunames),
                                                            dimnames = list(dmunames,dimnames(s$superslack_input.W)[[2]]))))
       supers.i.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(supers.i.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(s$superslack_input.W)[1])), 
                          data.frame(supers.i.Worst, row.names = NULL))
       colnames(supers.i.Worst)[3:(ncol(supers.i.Worst))] <- paste("superslack",
                                                         colnames(supers.i.Worst)[3:(ncol(supers.i.Worst))],
                                                         "Worst",sep = ".")
       
       supers.o.Worst <- do.call(rbind,lapply(seq(dim(s$superslack_output.W)[3]), 
                                         function(x) matrix(s$superslack_output.W[,,x],
                                                            nrow = length(dmunames),
                                                            dimnames = list(dmunames,dimnames(s$superslack_output.W)[[2]]))))
       supers.o.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(supers.o.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(s$superslack_output.W)[1])), 
                          data.frame(supers.o.Worst, row.names = NULL))
       colnames(supers.o.Worst)[3:(ncol(supers.o.Worst))] <- paste("superslack",
                                                         colnames(supers.o.Worst)[3:(ncol(supers.o.Worst))],
                                                         "Worst",sep = ".")
       
       supers.i.Best <- do.call(rbind,lapply(seq(dim(s$superslack_input.B)[3]), 
                                        function(x) matrix(s$superslack_input.B[,,x],
                                                           nrow = length(dmunames),
                                                           dimnames = list(dmunames,dimnames(s$superslack_input.B)[[2]])) ))
       supers.i.Best <- cbind(data.frame(DMU = dimnames(supers.i.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(s$superslack_input.B)[1])), 
                         data.frame(supers.i.Best, row.names = NULL))
       colnames(supers.i.Best)[3:(ncol(supers.i.Best))] <- paste("superslack",
                                                       colnames(supers.i.Best)[3:(ncol(supers.i.Best))],
                                                       "Best",sep = ".")
       supers.o.Best <- do.call(rbind,lapply(seq(dim(s$superslack_output.B)[3]), 
                                        function(x) matrix(s$superslack_output.B[,,x],
                                                           nrow = length(dmunames),
                                                           dimnames = list(dmunames,dimnames(s$superslack_output.B)[[2]])) ))
       supers.o.Best <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(supers.o.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(s$superslack_output.B)[1])), 
                         data.frame(supers.o.Best, row.names = NULL))
       colnames(supers.o.Best)[3:(ncol(supers.o.Best))] <- paste("superslack",
                                                       colnames(supers.o.Best)[3:(ncol(supers.o.Best))],
                                                       "Best",sep = ".")
       
       supers.i.df  <- merge(supers.i.Worst,supers.i.Best, by = c("alphacut","DMU"))
       # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
       srtidx <- (3:ncol(supers.i.df))
       srtidx <- t(matrix(srtidx, ncol = 2))
       dim(srtidx) <- c(1,length(srtidx))
       supers.i.df <- supers.i.df[,c(2,1,srtidx)]
       
       supers.o.df  <- merge(supers.o.Worst,supers.o.Best, by = c("alphacut","DMU"))
       # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
       srtidx <- (3:ncol(supers.o.df))
       srtidx <- t(matrix(srtidx, ncol = 2))
       dim(srtidx) <- c(1,length(srtidx))
       supers.o.df <- supers.o.df[,c(2,1,srtidx)]
       supers.df <- cbind(supers.i.df, supers.o.df[,3:ncol(supers.o.df)])
     }else {
       supers.df <- NULL
     }
   
     # Lambdas ------------------
     
     lmb <- lambdas(object)
     
     lmbmat.Worst <- do.call(rbind,lapply(seq(dim(lmb$Worst)[3]), function(x) lmb$Worst[,,x] ))
     lmbdf.Worst <- cbind(data.frame(DMU = dimnames(lmbmat.Worst)[[1]], 
                                     alphacut = rep(object$alpha,each = dim(lmb$Worst)[1])), 
                          data.frame(lmbmat.Worst, row.names = NULL))
     colnames(lmbdf.Worst)[3:(ncol(lmbdf.Worst))] <- paste("lambda",colnames(lmbdf.Worst)[3:(ncol(lmbdf.Worst))],"Worst",sep = ".")
     
     lmbmat.Best <- do.call(rbind,lapply(seq(dim(lmb$Best)[3]), function(x) lmb$Best[,,x] ))
     lmbdf.Best <- cbind(data.frame(DMU = dimnames(lmbmat.Best)[[1]], 
                                    alphacut = rep(object$alpha,each = dim(lmb$Best)[1])), 
                         data.frame(lmbmat.Best, row.names = NULL))
     colnames(lmbdf.Best)[3:(ncol(lmbdf.Best))] <- paste("lambda",colnames(lmbdf.Best)[3:(ncol(lmbdf.Best))],"Best",sep = ".")
     lmb.df  <- merge(lmbdf.Worst,lmbdf.Best, by = c("alphacut","DMU"))
     # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
     srtidx <- (3:ncol(lmb.df))
     srtidx <- t(matrix(srtidx, ncol = 2))
     dim(srtidx) <- c(1,length(srtidx))
     lmb.df <- lmb.df[,c(2,1,srtidx)]
     
     
     # Targets ---------------------
     tar <- targets(object)
     
     
     tar.i.Worst <- do.call(rbind,lapply(seq(dim(tar$target_input.W)[3]), 
                                       function(x) matrix(tar$target_input.W[,,x],
                                                          nrow = length(dmunames),
                                                          dimnames = list(dmunames,dimnames(tar$target_input.W)[[2]]))))
     tar.i.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(tar.i.Worst)[[1]], 
                                   alphacut = rep(object$alpha,each = dim(tar$target_input.W)[1])), 
                        data.frame(tar.i.Worst, row.names = NULL))
     colnames(tar.i.Worst)[3:(ncol(tar.i.Worst))] <- paste("target",
                                                       colnames(tar.i.Worst)[3:(ncol(tar.i.Worst))],
                                                       "Worst",sep = ".")
     
     tar.o.Worst <- do.call(rbind,lapply(seq(dim(tar$target_output.W)[3]), 
                                       function(x) matrix(tar$target_output.W[,,x],
                                                          nrow = length(dmunames),
                                                          dimnames = list(dmunames,dimnames(tar$target_output.W)[[2]]))))
     tar.o.Worst <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(tar.o.Worst)[[1]], 
                                   alphacut = rep(object$alpha,each = dim(tar$target_output.W)[1])), 
                        data.frame(tar.o.Worst, row.names = NULL))
     colnames(tar.o.Worst)[3:(ncol(tar.o.Worst))] <- paste("target",
                                                       colnames(tar.o.Worst)[3:(ncol(tar.o.Worst))],
                                                       "Worst",sep = ".")
     
     tar.i.Best <- do.call(rbind,lapply(seq(dim(tar$target_input.B)[3]), 
                                      function(x) matrix(tar$target_input.B[,,x],
                                                         nrow = length(dmunames),
                                                         dimnames = list(dmunames,dimnames(tar$target_input.B)[[2]])) ))
     tar.i.Best <- cbind(data.frame(DMU = dimnames(tar.i.Best)[[1]], 
                                  alphacut = rep(object$alpha,each = dim(tar$target_input.B)[1])), 
                       data.frame(tar.i.Best, row.names = NULL))
     colnames(tar.i.Best)[3:(ncol(tar.i.Best))] <- paste("target",
                                                     colnames(tar.i.Best)[3:(ncol(tar.i.Best))],
                                                     "Best",sep = ".")
     tar.o.Best <- do.call(rbind,lapply(seq(dim(tar$target_output.B)[3]), 
                                      function(x) matrix(tar$target_output.B[,,x],
                                                         nrow = length(dmunames),
                                                         dimnames = list(dmunames,dimnames(tar$target_output.B)[[2]])) ))
     tar.o.Best <- cbind(data.frame(DMU = object$data$dmunames,#dimnames(tar.o.Best)[[1]], 
                                  alphacut = rep(object$alpha,each = dim(tar$target_output.B)[1])), 
                       data.frame(tar.o.Best, row.names = NULL))
     colnames(tar.o.Best)[3:(ncol(tar.o.Best))] <- paste("target",
                                                     colnames(tar.o.Best)[3:(ncol(tar.o.Best))],
                                                     "Best",sep = ".")
     
     tar.i.df  <- merge(tar.i.Worst,tar.i.Best, by = c("alphacut","DMU"))
     # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
     srtidx <- (3:ncol(tar.i.df))
     srtidx <- t(matrix(srtidx, ncol = 2))
     dim(srtidx) <- c(1,length(srtidx))
     tar.i.df <- tar.i.df[,c(2,1,srtidx)]
     
     tar.o.df  <- merge(tar.o.Worst,tar.o.Best, by = c("alphacut","DMU"))
     # Re-arranging columns so they are sorted in col1.Worst col1.Best, col2.Worst col2.Best, .... 
     srtidx <- (3:ncol(tar.o.df))
     srtidx <- t(matrix(srtidx, ncol = 2))
     dim(srtidx) <- c(1,length(srtidx))
     tar.o.df <- tar.o.df[,c(2,1,srtidx)]
     tar.df <- cbind(tar.i.df, tar.o.df[,3:ncol(tar.o.df)])
     
     if(!modelkl %in% c("additive","addsupereff")){
     df <- cbind(eff.df, s.df[,3:ncol(s.df)], lmb.df[,3:ncol(lmb.df)], tar.df[,3:ncol(tar.df)])
     } else{
       df <- cbind(s.df[,3:ncol(s.df)], lmb.df[,3:ncol(lmb.df)], tar.df[,3:ncol(tar.df)])
     }
     if(modelkl %in% c("addsupereff","sbmsupereff")){
       df <- cbind(df,supers.df[3:ncol(supers.df)])
     }
     
     if(exportExcel){
       df.list <- list(efficiencies= eff.df, 
                       slacks = s.df,
                       superslacks = supers.df,
                       lambdas = lmb.df,
                       targets = tar.df)
       df.list[sapply(df.list, is.null)] <- NULL
       if(is.null(filename)){
         filename <- paste("ResultsDEA",Sys.time(),".xlsx", sep = "")
         filename <- gsub(" ","_",filename)
         filename <- gsub(":",".",filename)
       }
       write_xlsx(df.list, path = filename)
     }
     return(df)
 }


}