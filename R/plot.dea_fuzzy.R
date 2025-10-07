#' @title Plot for fuzzy DEA models.
#'   
#' @description Plot some attributes of a fuzzy DEA model (Guo-Tanaka, Kao-Liu  and
#' possibilistic models).
#' 
#' @param x An object of class \code{dea_fuzzy} obtained by a fuzzy DEA model function.
#' @param showPlots Logical. When TRUE (default) the plots are shown one by one. When it 
#' is FALSE the plots are not shown and are returned by the function (invisibly) as a 
#' list.
#' @param ... Ignored, for compatibility issues.
#' 
#' @returns Depending on the model, it returns some plots.
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
#' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking.
#' Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York.
#' \doi{10.1007/978-3-319-06647-9}
#'
#' @method plot dea_fuzzy
#' 
#' @importFrom igraph graph.adjacency degree "V<-" "V" plot.igraph
#' @importFrom ggplot2 ggplot geom_line geom_histogram geom_col scale_fill_gradientn
#' facet_wrap scale_x_discrete scale_y_continuous theme_bw theme scale_fill_identity
#' guides xlab ylab coord_flip aes ggtitle geom_bar geom_text geom_point geom_rect
#' scale_y_discrete
#' @importFrom methods show
#' @importFrom graphics plot
#' @importFrom stats runif complete.cases
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices hcl.colors
#' @import dplyr
#' 
#' @export

plot.dea_fuzzy <- function(x,
                           showPlots = TRUE, ...) {
  
  object <- x
  if (!is.dea_fuzzy(object)) {
    stop("Input should be of class dea_fuzzy!")
  }

 modelname <- object$modelname
 
 # For check CRAN pass...
 Period <- value <- DMUs <- Pos <- Count <- . <- index <- iseff <- ..count.. <- 
   eff.mean_eff <- Aspect <-  alpha_cut <- dL <- dR <- isefflab <- m <- variable <- Score <- 
   ID <- l <- u <- alpha <- `h-level` <- NULL
 `V<-` <- NULL
 
 eff <- efficiencies(object)
 
 if (modelname == "fuzzy_guotanaka") {
   
   EFF <- do.call(rbind, lapply(seq(dim(eff)[3]), function(x)
     eff[, , x]))
   DMUs <- rownames(EFF)
   rownames(EFF) <- NULL
   
   eff_df <- data.frame(cbind(data.frame(DMUs = DMUs, EFF)), row.names = NULL)
   
   eff_df$alpha_cut <- rep(paste("h = ", as.character(object$h)), each = length(object$dmu_eval))
   p <- eff_df %>% mutate(xmin = m - dL, xmax = m + dR) %>%
     pivot_longer(-c(DMUs, alpha_cut), names_to = "variable", values_to = "value") %>%
     filter(!variable %in% c("dL", "dR")) %>%
     ggplot(aes(x = value, y = DMUs, color = DMUs)) +
     geom_point(aes(shape = ifelse(variable == "m", "triangle", "circle"))) +
     geom_line() + facet_wrap( ~ alpha_cut) +
     theme_bw() + theme(legend.position = "none") +
     scale_y_discrete(limits = rev(rownames(eff))) +
     xlab("Score")
   
   show(p)
   
 } else if (grepl("kaoliu", modelname)) {

   eff_worst <- eff$Worst
   eff_best <- eff$Best
   label <- "Score"
   
   if (length(dim(eff_worst)) == 3) { # Models with several scores (nonradial and deaps, but not always)
     eff_worst <-
       do.call(cbind, lapply(seq(dim(eff$Worst)[3]), function(x)
         eff$Worst[, "mean_efficiency", x]))
     eff_best <-
       do.call(cbind, lapply(seq(dim(eff$Best)[3]), function(x)
         eff$Best[, "mean_efficiency", x]))
     label <- "Mean Score"
   }
   
   dmunames <- object$data$dmunames[object$dmu_eval]
   EFF <- data.frame(DMUs = rep(dmunames, each = length(object$alpha)))
   EFF$alpha <- rep(object$alpha, length(dmunames))
   EFF$l <- c(t(eff_worst))
   EFF$u <- c(t(eff_best))
   rango <- diff(range(c(eff_worst, eff_best)))
   p <- EFF %>% mutate(ID = rev(as.numeric(factor(DMUs)))) %>%
     ggplot() +
     geom_rect(aes(ymin = ID - 0.25,
                   ymax = ID + 0.25,
                   xmin = (l - rango / 1000),
                   xmax = (u + rango / 1000),
                   fill = alpha)) +
     scale_y_continuous(breaks = seq_along(unique(EFF$DMUs)),
                        labels = rev(unique(EFF$DMUs)))+
     scale_fill_gradientn(colours = hcl.colors(10, palette = "YlGnBu", rev = TRUE)) +
     xlab(label) + ylab("DMUs") + 
     theme_bw()
   
   show(p)
   
 } else { # Possibilistic model
   
   p <- eff %>% as.data.frame() %>% cbind(data.frame(DMUs = rownames(.)),.) %>% 
     pivot_longer(-1, names_to = "h-level", values_to = "Score") %>% 
     mutate(`h-level` = as.numeric(`h-level`)) %>%
     ggplot(aes(x = `h-level`, y = Score, colour = DMUs)) +  
     geom_line() + geom_point() + theme_bw()
   
    show(p)
    
 }
 
 invisible(list(`Efficiencies` = p))

}
