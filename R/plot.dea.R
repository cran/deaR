#' @title Plot for DEA models.
#'   
#'   
#' @description Plot some attribute of a DEA model (conventional, fuzzy or Malmquist).
#' 
#' @param x An object of class \code{"dea"} obtained by a dea model function.
#' @param showPlots Logical. When TRUE (default) the plots are shown one by one. When it 
#' is FALSE the plots are not shown and are returned by the function (invisiblily) as a 
#' list.
#' @param ... Ignored, for compatibility issues.
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
#' data_example <- read_data(datadea = Fortune500,
#'                           inputs = 2:4, 
#'                           outputs = 5:6)
#' result <- model_basic(data_example)
#' plot(result)
#' 
#' @references
#' #' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking. Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York. DOI: 10.1007/978-3-319-06647-9
#' @method plot dea
#' @importFrom igraph graph.adjacency degree "V<-" "V" plot.igraph
#' @importFrom ggplot2 ggplot geom_line geom_histogram geom_col facet_wrap scale_x_discrete theme_bw scale_fill_identity guides xlab ylab coord_flip aes ggtitle geom_bar geom_text
#' @importFrom methods show
#' @importFrom graphics plot
#' @importFrom grDevices pdf
#' @importFrom stats runif complete.cases
#' @importFrom gridExtra grid.arrange
#' @import plotly dplyr
#' @export
#' 

plot.dea <- function(x, showPlots = TRUE, ...){
  object <- x
  orientation <- x$orientation
  if (!is.dea(object)) {
    stop("Input should be of class dea!")
  }

 modelname <- object$modelname
 # For check CRAN pass...
 Period <- value <- DMU <- Pos <- Count <- . <- index <- iseff <- isefflab <- ..count.. <- eff.mean_eff <- Aspect <-  NULL
 `V<-` <- NULL
 if (modelname == "deaps") {
   if (object$restricted_eff == FALSE) {
     warning("Plotting a Preference Structure model with unrestricted efficiencies is not
             available!")
     return(NULL)
   }
 }
 if (grepl("add", modelname)) {
   warning("Plotting additive models are not implemented yet!")
   return(NULL)
 }
 if (modelname == "basic") {
   if (!object$orientation %in% c("io", "oo")) {
     warning("Plotting Basic model with is only available with input/output orientations!")
     return(NULL)
   }
 } 
 
 
 if(modelname %in% c("malmquist", "cross_efficiency")){
   if (modelname == "malmquist") {
     malmdata <- summary(object, exportExcel = FALSE)
     results <- malmdata$Results
     sumres <- malmdata$means_by_period
     #colnames(results) <- c("Period", "DMU", "EC", "TC", "PTEC", "SEC", "MI")
     #colnames(sumres) <- c("Period", "EC", "TC", "PTEC", "SC", "MI")
     #colnames(results) <- toupper(colnames(results))
     #colnames(sumres) <- toupper(colnames(sumres))
     #results <- results[,c("Period", "DMU", "EC", "PTEC", "SEC", "TC", "MI")]
     
     
     
     results %>% gather(key = "index", value = "value", -c("Period", "DMU")) -> resmelt 
     resmelt$index <- as.factor(resmelt$index)
     resmelt$index <- factor(resmelt$index, levels(resmelt$index)[c(1,3,4,5,2)])
     resmelt %>% group_by(index) %>% 
       do(p = plot_ly(., x =~Period, y=~value, color = ~DMU, mode = 'lines', type = 'scatter', colors = "Paired") %>% 
            layout(yaxis = list(title=~index)))  %>% 
       subplot(nrows = NROW(.), shareX = TRUE,titleY=TRUE)-> resplot
     
     
     
     sumres %>% gather(key= "index", value = "value", -c("Period")) -> sumresmelt 
     sumresmelt$index <- as.factor(sumresmelt$index)
     sumresmelt$index <- factor(sumresmelt$index, levels(sumresmelt$index)[c(1,3,4,5,2)])
     sumresmelt %>% plot_ly(x = ~Period, y = ~value, type = 'scatter', mode = 'lines', color = ~index) -> sumplot
     
     
     if (showPlots) {
        invisible(readline(prompt = "Press [enter] to continue"))
        pdf(NULL)
        show(resplot)
        invisible(readline(prompt = "Press [enter] to continue"))
        pdf(NULL)
        show(sumplot)
     }
     invisible(list(`Results plot` = resplot, `Summary plot` = sumplot))
   } else {
      results <- list(Arb = object$Arbitrary$cross_eff,
                      M2_agg = object$M2_agg$cross_eff,
                      M2_ben = object$M2_ben$cross_eff,
                      M3_agg = object$M3_agg$cross_eff,
                      M3_ben = object$M3_ben$cross_eff)
      titles <- list("Arbitrary Method", "Method II - Aggresive", "Method II - Benevolent",
                     "Method III - Aggresive", "Method III - Benevolent")
      titles[sapply(results, is.null)] <- NULL 
      results[sapply(results, is.null)] <- NULL
      p <- list()
      for (i in seq_along(results)) {
         xlab <- ylab <- colnames(results[[i]])
         p[[i]] <- plot_ly(x = xlab, y = rev(ylab), z  = results[[i]][nrow(results[[i]]):1, ], colors = "RdYlGn",type = "heatmap" ) %>% 
            layout(title = titles[[i]])
         if (showPlots) {
            invisible(readline(prompt = "Press [enter] for next plot"))
            pdf(NULL)
            show(p[[i]])
         }
      }
      names(p) <- titles
      invisible(p)
      
   }
 }else{
   
   # Efficiencies histogram ----------
   if (!modelname %in% c("nonradial", "deaps")) {
     eff <-
       data.frame(DMU = object$data$dmunames[object$dmu_eval], eff = efficiencies(object))
     eff <- eff[complete.cases(eff), ]
     
     if (modelname == "profit") {
       # Efficient DMUS have eff = 1 and optimio  = io 
       iseff <- sapply(seq_len(length(object$DMU)), function(i) {
         sum(c((object$DMU[[i]]$optimal_output - object$data$output[, i]) ^ 2 , 
               (object$DMU[[i]]$optimal_input - object$data$input[, i]) ^ 2)) < 1e-10
       }
       )
       eff$iseff <- ifelse(iseff, 1, 0)
       
       
     } else {
     # Efficient DMUS have eff = 1 and slacks  = 0
     slk <- slacks(object)
     slk[sapply(slk, is.null)] <- NULL
     null_slk <-
       apply(
         do.call(cbind, slk),
         MARGIN = 1,
         FUN = function(x)
           sum(x ^ 2)
       ) < 1e-4 # Null slacks
     eff_1 <- abs(eff$eff - 1) < 1e-8 # Efficiency = 1
     
     eff$iseff <-
        ifelse(null_slk & eff_1, # Slacks = 0 and efficiency = 1
               1, 
               ifelse(eff$eff >= 1 + 1e-8 & orientation != "oo", 
                      2, 
                      0)
               )
     
     
     eff$isefflab <- ifelse(eff$iseff == 1, "Efficient",ifelse(eff$iseff == 0, "Inefficient","Super-efficient"))
     }
     if (!modelname %in% c("supereff_basic", "sbmsupereff")) {
        
        
        if (sum(eff$iseff != 1)  > 0){
           eff %>% filter(iseff != 1) %>% ggplot(aes(x = eff)) +
              geom_histogram(breaks = c(seq(
                 from = min(eff$eff),
                 to = max(eff$eff[eff$iseff != 1]),
                 length.out = 10
              )),
              col = "white",
              aes(fill = ifelse(iseff == 0, "red", ifelse(iseff == 1,"lightgreen","lightblue")))) +
              theme_bw() + scale_fill_identity() + xlab("Efficiency") + ylab("Count") +
              ggtitle("Non-efficient DMU distribution") -> p2
        } else{
           warning("All DMUs are efficient!")
        }
       
       eff %>% ggplot(aes(x = isefflab)) + 
          geom_bar(aes(fill = ifelse(iseff == 0, "red", ifelse(iseff == 1,"lightgreen","lightblue")))) +
          theme_bw() + scale_fill_identity() + xlab("") + #scale_x_discrete(labels = isefflab) +
          ylab("Count") +
         ggtitle("Efficient/Non Efficient DMUs") + geom_text(stat = 'count',
                                                             aes(label = ..count..),
                                                             vjust = -0.5) -> p1
       if (showPlots) {
          #pdf(NULL)
          if(sum(eff$iseff != 1)  > 0){
         grid.arrange(p1, p2, nrow = 1)
          } else{
             grid.arrange(p1, nrow = 1)
          }
       }
     } else {
        eff %>% mutate(isefflabel = ifelse(iseff == 1, "Efficient",ifelse(iseff == 0, "Inefficient","Superefficient"))) %>% 
           filter(iseff !=1) %>% ggplot(aes(x = eff)) + 
           geom_histogram(bins = 5,
                          col = "white",
                          aes(fill = ifelse(iseff == 0,
                                            "red",
                                            ifelse(iseff == 1,
                                                   "lightgreen","lightblue")))) +
           theme_bw() + scale_fill_identity() + xlab("Efficiency / Superefficiency") + ylab("Count") +
           facet_wrap(~isefflabel, scales = "free")-> p2
           p1 <- NULL
           
        if (showPlots) {
          # pdf(NULL)
           show(p2)
        }
     }
     
   
   } else {
     eff <-
       data.frame(DMU = object$data$dmunames, eff = efficiencies(object))
     eff %>% mutate(iseff = ifelse(abs(eff.mean_eff - 1) < 1e-4, 1 , 0)) -> eff
     effmelted <-
       eff %>% gather(key = "Aspect", value = "eff",-c(DMU, iseff))
     effmelted$Aspect <- gsub("eff.", "", x = effmelted$Aspect)
     effmelted %>% filter(iseff == 0 & Aspect != "mean_eff") %>%
       ggplot(aes(x = eff)) + geom_histogram(bins = 10, fill = "red", col = "white") +
       theme_bw() + scale_fill_identity() + xlab("Efficiency") + ylab("Count") + facet_wrap( ~
                                                                                               Aspect, scales = "free") -> p2
     
     
     effmelted %>% filter(Aspect == "mean_eff") %>% ggplot(aes(x = factor(iseff))) +
       geom_bar(aes(fill = ifelse(iseff == 0, "red", "lightgreen")), width = 0.5) + theme_bw() +
       scale_fill_identity() + xlab("") + scale_x_discrete(labels = c("Non-efficient", "Efficient")) + ylab("Count") +
       ggtitle("Efficient/Non Efficient DMUs") + geom_text(stat = 'count', aes(label =
                                                                                 ..count..), vjust = -0.5) -> p1
     
     if (showPlots) {
    # pdf(NULL)
     show(p1)
     invisible(readline(prompt = "Press [enter] to continue"))
    # pdf(NULL)
     show(p2)
     }
   }
   
   # Reference ranking --------------
 
   if ((!modelname %in% c("supereff_basic", "sbmsupereff")) & (sum(eff$iseff != 1)  > 0)) {
     ref <- references(object)
     lmbd <- lambdas(object)
     dmunames <- object$data$dmunames[object$dmu_eval]
     refnames <- unique(unlist(lapply(ref, function (x)
       names(x))))
     urefnames <- names(ref)
     effdmus <- dmunames[which(!dmunames %in% urefnames)]
     refnames <- refnames[refnames %in% effdmus]
     
     RefMat <-
       matrix(
         0,
         nrow = length(urefnames),
         ncol = length(effdmus),
         dimnames = list(urefnames, sort(effdmus))
       )
     RefMat[urefnames, refnames] <- round(lmbd[urefnames, refnames], 4)
     RefMatl <- RefMat > 0
     ranking <- sort(colSums(RefMatl))
     realcount <- ranking
     ranking[ranking == 0] <- min(ranking[ranking > 0]) * 1e-2
     effnames <- names(ranking)
     
     rkdf <-
       data.frame(
         DMU = effnames,
         Count = ranking,
         realcount = realcount,
         Pos = seq(length(effdmus))
       )
     
     rkdf %>% ggplot(aes(x = as.factor(Pos), y = Count)) + geom_col(aes(fill = Count)) + coord_flip() +
       scale_x_discrete(labels = rkdf$DMU, drop = FALSE) + theme_bw() + xlab("Efficient DMUs") +
       ylab("# times appearing in reference sets") + geom_text(aes(label = realcount), hjust =
                                                                 -0.5) +
       guides(fill = FALSE) -> refplot
     
     if (showPlots){
     #  pdf(NULL)
       invisible(readline(prompt = "Press [enter] to continue"))
       show(refplot)
       invisible(readline(prompt = "Press [enter] to continue"))
     }
   } else {
      refplot <- NULL
     warning("Ranking plots with those models are not implemented yet!")
     
     
   }
    
   
   # Reference Graph ------------------------
   
   if ((!modelname %in% c("supereff_basic", "sbmsupereff"))&(sum(eff$iseff != 1)  > 0)) {
      #dimlmbd <- dim(lmbd)
     lmbd <- matrix(lmbd[complete.cases(lmbd), ], ncol = ncol(lmbd), dimnames = dimnames(lmbd))
    
     
     # Make the matrix square in case dmu_eval != dmu_ref
     
     allnames <- unique(unlist(dimnames(lmbd)))
     
     extendedlmbd <- matrix(0,
                            nrow = length(allnames),
                            ncol = length(allnames),
                            dimnames = list(allnames, allnames))
     extendedlmbd[dimnames(lmbd)[[1]],dimnames(lmbd)[[2]]] <- lmbd
     
     adjmatrix <- extendedlmbd > 0
     G <- graph.adjacency(adjmatrix, diag = FALSE)
     
     efficient <- which(dmunames %in% effdmus)
     non_efficient <- which(dmunames %in% urefnames) # Non-efficient includes inefficent and superefficient DMUs
     names(non_efficient) <- ifelse(eff$iseff[non_efficient] == 0,"red","lightblue")
     
     rnd <- runif(1)
     nefflocX <-
       2 * cos(2 * pi * (0:(length(non_efficient) - 1)) / length(non_efficient) + rnd)
     nefflocY <-
       2 * sin(2 * pi * (0:(length(non_efficient) - 1)) / length(non_efficient) + rnd)
     efflocX <-
       4 * cos(2 * pi * (0:(length(efficient) - 1)) / length(efficient))
     efflocY <-
       4 * sin(2 * pi * (0:(length(efficient) - 1)) / length(efficient))
     
     relations <- degree(G, mode = "in") - degree(G, mode = "out")
     
     
     locX <- numeric(length(relations))
     locY <- numeric(length(relations))
     locX[efficient] <- efflocX
     locX[non_efficient] <- nefflocX
     
     locY[efficient] <- efflocY
     locY[non_efficient] <- nefflocY
     locations <- cbind(locX , locY)
     
     V(G)$color <- numeric(length(dmunames))
     V(G)$color[efficient] <- rep("green", length(efficient))
     V(G)$color[non_efficient] <- names(non_efficient) #rep("red", length(non_efficient))
     lmbd2 <- lmbd / rowSums(lmbd)
     lmbd2[is.nan(lmbd2)] <- 0
     V(G)$size <-
       colSums(lmbd2) ^ 1.1 + 10#0.6*degree(G, mode = "in") + 10
     
     
     graphPlot <- list(G = G, locations = locations)  
    
     if (showPlots) {
      # pdf(NULL)
       plot(
         G,
         layout = locations,
         xlim = c(-2, 2),
         ylim = c(-.4, .4),
         edge.arrow.size = 0.75,
         edge.arrow.width = 0.25,
         edge.curved = FALSE
       )
     }
   }  else {
   graphPlot <- NULL
   warning("Graph plots with those models are not implemented yet!")
   
   
}
   invisible(list(`Eff/Ineff count` = p1, `Ineff Dstr` = p2, `References Plot`= refplot, `References Graph` = graphPlot))
 }
 
 
}
