#' @title Plot for DEA models.
#'   
#' @description Plot some attribute of a DEA model.
#' 
#' @param x An object of class \code{dea} obtained by a DEA model function.
#' @param tol Numeric. Absolute tolerance for numeric comparisons. By default, it is 1e-4.
#' @param showPlots Logical. When TRUE (default) the plots are shown one by one. When it 
#' is FALSE the plots are not shown and are returned by the function (invisibly) as a list.
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
#' @examples
#' data_example <- make_deadata(datadea = Fortune500,
#'                              inputs = 2:4, 
#'                              outputs = 5:6)
#' result <- model_basic(data_example)
#' plot(result)
#' 
#' @references
#' Zhu, J. (2014). Quantitative Models for Performance Evaluation and Benchmarking.
#' Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer, New York.
#' \doi{10.1007/978-3-319-06647-9}
#'
#' @method plot dea
#'
#' @importFrom igraph graph.adjacency degree "V<-" "V" plot.igraph
#' @importFrom ggplot2 ggplot geom_col theme_bw geom_histogram scale_x_discrete
#' scale_fill_identity guides xlab ylab coord_flip aes ggtitle geom_text theme
#' scale_fill_manual
#' @importFrom methods show
#' @importFrom graphics plot
#' @importFrom grDevices pdf
#' @importFrom stats runif complete.cases
#' @importFrom gridExtra grid.arrange
#' @import plotly dplyr
#' 
#' @export

plot.dea <- function(x,
                     tol = 1e-4,
                     showPlots = TRUE, ...) {
  
  deasol <- x
  
  if (!is.dea(deasol)) {
    stop("Input should be of class dea. Run a model first!")
  }
  
  modelname <- deasol$modelname
  
  p1 <- p2 <- refplot <- graphplot<- NULL
  
  if (!modelname %in% c("additive", "addsupereff", "basic", "deaps", "fdh_basic", "multiplier",
                        "nonradial", "profit", "rdm", "sbmeff", "sbmsupereff", "supereff_basic",
                        "cross_efficiency")) {
    stop("Model not supported.")
  }
  
  # For check CRAN pass...
  Score <- group <- NULL  
  
  
  if (modelname == "cross_efficiency") {
    # Cross Efficiency -----------------
    results <- list(Arb = deasol$Arbitrary$cross_eff,
                    M2_agg = deasol$M2_agg$cross_eff,
                    M2_ben = deasol$M2_ben$cross_eff,
                    M3_agg = deasol$M3_agg$cross_eff,
                    M3_ben = deasol$M3_ben$cross_eff)
    titles <- list("Arbitrary Method", "Method II - Aggressive", "Method II - Benevolent",
                   "Method III - Aggressive", "Method III - Benevolent")
    titles[sapply(results, is.null)] <- NULL 
    results[sapply(results, is.null)] <- NULL
    p <- list()
    for (i in seq_along(results)) {
      xlab <- ylab <- colnames(results[[i]])
      p[[i]] <- plot_ly(x = xlab, y = rev(ylab), z  = results[[i]][nrow(results[[i]]):1, ],
                        colors = "RdYlGn",type = "heatmap" ) %>% 
        layout(title = titles[[i]])
      if (showPlots) {
        if (i > 1) invisible(readline(prompt = "Press [enter] to continue"))
        #pdf(NULL)
        print(p[[i]])
      }
    }
    names(p) <- titles
    invisible(p)
        
    # } else if (modelname == "malmquist") {
    # #   Malmquist ------------------------
    #
    #   malmdata <- summary(deasol, exportExcel = FALSE)
    #   results <- malmdata$Results
    #   sumres <- malmdata$means_by_period
    # 
    #   results %>% gather(key = "index", value = "value", -c("Period", "DMU")) -> resmelt
    #   resmelt$index <- as.factor(resmelt$index)
    #   resmelt$index <- factor(resmelt$index, levels(resmelt$index)[c(1, 3, 4, 5, 2)])
    #   resmelt %>% group_by(index) %>%
    #     do(p = plot_ly(., x = ~Period, y = ~value, color = ~DMU, mode = 'lines',
    #                    type = 'scatter', colors = "Paired") %>%
    #          layout(yaxis = list(title = ~index))) %>%
    #     subplot(nrows = NROW(.), shareX = TRUE, titleY = TRUE) -> resplot
    # 
    #   # sumres %>% gather(key = "index", value = "value", -c("Period")) -> sumresmelt
    #   # sumresmelt$index <- as.factor(sumresmelt$index)
    #   # sumresmelt$index <- factor(sumresmelt$index, levels(sumresmelt$index)[c(1, 3, 4, 5, 2)])
    #   # sumresmelt %>% plot_ly(x = ~Period, y = ~value, type = 'scatter', mode = 'lines',
    #   #                        color = ~index) -> sumplot
    # 
    #   if (showPlots) {
    #     pdf(NULL)
    #     show(resplot)
    #     # invisible(readline(prompt = "Press [enter] to continue"))
    #     # pdf(NULL)
    #     # show(sumplot)
    #   }
    #   invisible(resplot)
    #   #invisible(list(`Results plot` = resplot, `Summary plot` = sumplot))
    
  } else {
    
    ######## Scores histogram ######## -------------------
    
    if (modelname %in% c("addsupereff", "sbmsupereff", "supereff_basic")) {
      super <- TRUE
    } else {
      super <- FALSE
    }
    
    if (modelname %in% c("nonradial", "deaps"))
    {
      xlabel <- "Mean score"
    } else {
      xlabel <- "Score"
    }
    
    dmu_eval <- deasol$dmu_eval
    
    # Where is the efficiency score in the deasol$DMU[[1]] list?
    if (modelname == "profit") {
      iesc <- 2
    } else {
      iesc <- 1
    }
    eff <- unlist(lapply(deasol$DMU, function(x) x[[iesc]]))
    i_nona <- which(!is.na(eff)) # Indices in dmu_eval which are not NA in eff
    
    dmu_eff <- eff_dmus(deasol, tol = tol) # Efficient DMUs
    i_eff <- which(dmu_eval %in% dmu_eff) # Indices in dmu_eval of efficient DMUs
    ndeff <- length(i_eff)
    i_ineff <- i_nona[which(!i_nona %in% i_eff)] # Indices in dmu_eval of inefficient DMUs
    ndineff <- length(i_ineff)
    
    if (modelname == "supereff_basic") {
      # Inefficient histogram
      if (ndineff > 0) {
        df <- data.frame(Score = eff[i_ineff])
        p1 <- ggplot(df, aes(x = Score)) + 
          geom_histogram(breaks = c(seq(from = min(eff[i_ineff]) - tol,
                                        to = max(eff[i_ineff]) + tol,
                                        length.out = 10)),
                         col = "white", fill = "red") + 
          theme_bw() + scale_fill_identity() + xlab(xlabel) + ylab("Count") +
          ggtitle("Inefficient DMUs distribution")
      } 
    } else {
      # Number of efficient/inefficient histogram
      df <- data.frame(group = c("Efficient", "Inefficient"),
                       count = c(ndeff, ndineff))
      p1 <- ggplot(df, aes(x = group, y = count, fill = group)) + 
        geom_col() + theme_bw() +
        scale_fill_manual(values = c("Efficient" = "lightgreen", "Inefficient" = "red")) +
        xlab(" ") + ylab("Count") + ggtitle("Efficient/Inefficient DMUs") +
        geom_text(aes(label = count), vjust = -0.5) +
        theme(legend.position = "none")
    }
    
    if (super) {
      # Efficient histogram
      df <- data.frame(Score = eff[i_eff])
      if (ndeff > 0) {
        p2 <- ggplot(df, aes(x = Score)) + 
          geom_histogram(breaks = c(seq(from = min(eff[i_eff]) - tol,
                                        to = max(eff[i_eff]) + tol,
                                        length.out = 10)),
                         col = "white", fill = "lightgreen") + 
          theme_bw() + scale_fill_identity() + xlab(xlabel) + ylab("Count") +
          ggtitle("Efficient DMUs distribution")
      }
    }else {
      # Inefficient histogram
      if (ndineff > 0) {
        df <- data.frame(Score = eff[i_ineff])
        p2 <- ggplot(df, aes(x = Score)) + 
          geom_histogram(breaks = c(seq(from = min(eff[i_ineff]) - tol,
                                        to = max(eff[i_ineff]) + tol,
                                        length.out = 10)),
                         col = "white", fill = "red") + 
          theme_bw() + scale_fill_identity() + xlab(xlabel) + ylab("Count") +
          ggtitle("Inefficient DMUs distribution")
      }
    }
    
    if (showPlots) {
      if(is.null(p1)){
        grid.arrange(p2, nrow = 1)
      }else if(is.null(p2)){
        grid.arrange(p1, nrow = 1)
      }else{
        grid.arrange(p1, p2, nrow = 1)
      }
    }
    
    ######## References plots ######## ----------------
    
    if ((!super) && (ndineff != 0)) {
      
      ##### References Ranking ##### ----------
      
      dmu_ref <- deasol$dmu_ref
      ndr <- length(dmu_ref)
      lmbnum <- matrix(lambdas(deasol)[i_ineff, ], ncol = ndr)
      lmb <- (lmbnum > tol)
      colnames(lmb) <- names(dmu_ref)
      rownames(lmb) <- names(i_ineff)
      
      nref <- dmu_ref
      for (i in 1:ndr) {
        nref[i] <- sum(lmb[, i])
      }
      i_ref <- which(nref > 0) # Indices in dmu_ref of reference DMUs (references of inefficient DMUs)
      aux <- which(dmu_ref %in% dmu_eff)
      i_effnoref <- aux[which(!aux %in% i_ref)] # Indices in dmu_ref of efficient DMUs that are not reference
      
      datos <- sort(nref[c(i_ref, i_effnoref)])
      df <- data.frame(group = names(datos),
                       count = datos)
      refplot <- ggplot(df, aes(x = group, y = count, fill = count)) + 
        geom_col() + theme_bw() + coord_flip() +
        scale_x_discrete(limits = df$group, drop = FALSE) +
        xlab("Efficient DMUs") + ylab("# times appearing in reference sets") +
        geom_text(aes(label = count), hjust = -0.5) +
        theme(legend.position = "none")
      
      if (showPlots) {
        invisible(readline(prompt = "Press [enter] to continue"))
        show(refplot)
      }
      
      ##### References Graph ##### ---------------
      
      ndref <- length(i_ref)
      ndeffnoref <- length(i_effnoref)
      ndeff <- ndref + ndeffnoref
      nda <- ndineff + ndeff
      names_all <- c(names(dmu_eval[i_ineff]), names(dmu_ref[c(i_ref, i_effnoref)]))
      
      # Adjacency matrix
      Adj <- matrix(FALSE, nrow = nda, ncol = nda, dimnames = list(names_all, names_all)) 
      Adjnum <- Adj
      if ((ndineff > 0) && (ndref > 0)) {
        Adj[1:ndineff, (ndineff + 1):(ndineff + ndref)] <- lmb[, i_ref]
        Adjnum[1:ndineff, (ndineff + 1):(ndineff + ndref)] <-
          (lmbnum[, i_ref] / rowSums(matrix(lmbnum[, i_ref], nrow = ndineff)))
      }
      
      rnd <- 0.12345678
      inefflocX <- 2 * cos(2 * pi * (0:(ndineff - 1)) / ndineff + rnd)
      inefflocY <- 2 * sin(2 * pi * (0:(ndineff - 1)) / ndineff + rnd)
      efflocX <- 4 * cos(2 * pi * (0:(ndeff - 1)) / ndeff)
      efflocY <- 4 * sin(2 * pi * (0:(ndeff - 1)) / ndeff)
      locations <- cbind(c(inefflocX, efflocX), c(inefflocY, efflocY))
      
      G <- graph.adjacency(Adj)
      V(G)$color <- c(rep("red", ndineff), rep("lightgreen", ndeff))
      V(G)$size <- colSums(Adjnum) ^ 1.1 + 10
      
      graphplot <- list(G = G, locations = locations)
      
      if (showPlots) {
        invisible(readline(prompt = "Press [enter] to continue"))
        plot(G, layout = locations, xlim = c(-2, 2), ylim = c(-.4, .4),
             edge.arrow.size = 0.75, edge.arrow.width = 0.25, edge.curved = FALSE)
      }
    }
    invisible(list(`Graph 1` = p1, `Graph 2` = p2, `References Ranking`= refplot,
                   `References Graph` = graphplot))
  }
  
}