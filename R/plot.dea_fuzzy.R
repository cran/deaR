#' @title Plot for Fuzzy DEA models.
#'   
#'   
#' @description Plot some attribute of a Fuzzy DEA model (Guo-Tanaka and Kao-Liu models).
#' 
#' @param x An object of class \code{"dea_fuzzy"} obtained by a fuzzy dea model function.
#' @param showPlots Logical. When TRUE (default) the plots are shown one by one. When it 
#' is FALSE the plots are not shown and are returned by the function (invisiblily) as a 
#' list.
#' @param ... Ignored, for compatibility issues.
#' 
#'   
#' @return Depending on the model it returns ...       
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
#' 
#' @references
#' #' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking. Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York. DOI: 10.1007/978-3-319-06647-9
#' @method plot dea_fuzzy
#' @importFrom igraph graph.adjacency degree "V<-" "V" plot.igraph
#' @importFrom ggplot2 ggplot geom_line geom_histogram geom_col scale_fill_gradientn facet_wrap scale_x_discrete scale_y_continuous theme_bw scale_fill_identity guides xlab ylab coord_flip aes ggtitle geom_bar geom_text geom_point geom_rect
#' @importFrom methods show
#' @importFrom graphics plot
#' @importFrom stats runif complete.cases
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices hcl.colors
#' @import dplyr
#' @export
#' 

plot.dea_fuzzy <- function(x, showPlots = TRUE, ...) {
  object <- x
  if (!is.dea_fuzzy(object)) {
    stop("Input should be of class dea_fuzzy!")
  }

 modelname <- object$modelname
 # For check CRAN pass...
 Period <- value <- DMU <- Pos <- Count <- . <- index <- iseff <- ..count.. <- 
   eff.mean_eff <- Aspect <-  alpha_cut <- dL <- dR <- isefflab <- m <- variable <- Efficiency <- 
   ID <- l <- u <- alpha <- `h-value` <- NULL
 `V<-` <- NULL
 
 
 # Checking orientation
 if(!object$orientation %in% c("io","oo")){
   stop("Orientation should be either io or oo!")
 }
 
 eff <- efficiencies(object)
 # If output oriented we obtain efficiencies from scores
 if(object$orientation == "oo"){
   eff <- lapply(eff, function(x) 1/x)
 }
 
 if (modelname == "fuzzy_guotanaka") {
   EFF <- do.call(rbind, lapply(seq(dim(eff)[3]), function(x)
     eff[, , x]))
   DMU <- dimnames(EFF)[[1]]
   rownames(EFF) <- NULL
   
   eff_df <-
     data.frame(cbind(data.frame(DMU = DMU, EFF)), row.names = NULL)
   
   eff_df$alpha_cut <- rep(as.numeric(dimnames(eff)[[3]]), each = length(object$dmu_eval))
   p <- eff_df %>% mutate(xmin = m - dL, xmax = m + dR) %>%
     gather(variable, value,-c(DMU, alpha_cut)) %>%
     filter(!variable %in% c("dL", "dR")) %>%
     ggplot(aes(x = value, y = DMU, color = DMU)) + geom_point() +
     geom_line() +
     facet_wrap( ~ alpha_cut) +
     theme_bw() +
     xlab("Efficiency")
   show(p)
   
 } else if (grepl("kaoliu",modelname)) {
   if (grepl("add", modelname)) {
     warning("Plotting additive models are not implemented yet!")
     return(NULL)
   }
   # Checking for number of efficiencies (radial vs nonradial models)
   neff <- length(object$alphacut[[1]]$DMU$Worst[[1]]$efficiency) 
   
   
   if (neff == 1){ # radial model
   eff_worst <- eff$Worst
   eff_best <- eff$Best
   dmunames <- object$data$dmunames
   cnames <- unlist(dimnames(eff_worst)[2])
   L <- eff$Worst
   U <- eff$Best
   label <- "Efficiency"
   } else { # non-radial model
     cnames <- dimnames(eff$Worst)[[3]]
     eff_worst <-
       do.call(cbind, lapply(seq(dim(eff$Worst)[3]), function(x)
         eff$Worst[,"mean_efficiency", x]))
     eff_best <-
       do.call(cbind, lapply(seq(dim(eff$Best)[3]), function(x)
         eff$Best[,"mean_efficiency", x]))
     L <- eff_worst
     U <- eff_best
     label <- "Mean Efficiency"
   }
   
   
   
   EFF <- data.frame(DMU = rep(object$data$dmunames, each= length(object$alpha)))
   EFF$alpha <- rep(object$alpha, length(object$data$dmunames))
   EFF$l <- c(t(L))
   EFF$u <- c(t(U))
   p <- EFF %>% mutate(ID = rev(as.numeric(factor(DMU)))) %>%
     ggplot() +
     geom_rect(aes(ymin = ID - (1/4),
                   ymax = ID + (1/4),
                   xmin = (l-l/1000),
                   xmax = (u+u/1000),
                   fill = alpha)) +
     scale_y_continuous(breaks = seq_along(unique(EFF$DMU)),
                        labels = rev(unique(EFF$DMU)))+
     scale_fill_gradientn(colours = hcl.colors(10, palette = "YlGnBu", rev = TRUE)) +
     xlab(label) + ylab("DMU") + 
     theme_bw()
   show(p)
   
   
   
   
   # rownames(eff_worst) <- NULL
   # rownames(eff_best) <- NULL
   # eff_worst <-
   #   data.frame(DMU = dmunames, eff_worst, row.names = NULL)
   # eff_best <-
   #   data.frame(DMU = dmunames, eff_best, row.names = NULL)
   # colnames(eff_worst) <- c("DMU", cnames)
   # colnames(eff_best) <- c("DMU", cnames)
   # eff_worst <- gather(eff_worst, variable, value,-c(DMU))
   # eff_best <- gather(eff_best, variable, value, -c(DMU))
   # colnames(eff_worst) <- c("DMU", "alpha_cut", ifelse(neff == 1, "Efficiency", "Mean Efficiency"))
   # colnames(eff_best) <- c("DMU", "alpha_cut", ifelse(neff == 1, "Efficiency", "Mean Efficiency"))
   # eff_worst$Case <- "Worst"
   # eff_best$Case <- "Best"
   # eff_df <- rbind(eff_worst, eff_best)
   # eff_df$DMU <- as.factor(eff_df$DMU)
   # p <- eff_df  %>% ggplot(aes(
   #   x = Efficiency,
   #   y = DMU,
   #   color = DMU)) +
   #   geom_point() + geom_line() + ylab("DMU") + 
   #   facet_wrap( ~ alpha_cut)
 } else {
   p <- eff %>% as.data.frame() %>% cbind(data.frame(DMU = rownames(.)),.) %>% 
     pivot_longer(-1, names_to = "h-value", values_to = "Efficiency") %>% 
     mutate(`h-value` = as.numeric(`h-value`)) %>%
     ggplot(aes(x = `h-value`, y = Efficiency, colour = DMU)) +  
     geom_line() + geom_point() + theme_bw()
    show(p)
 }
 # if (showPlots) {
 #   ggplotly(p)
 # }else{
 #   show(p)
 # }
 
 invisible(list(`Efficiencies` = p))
}