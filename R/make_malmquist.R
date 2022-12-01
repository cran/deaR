#' @title make_malmquist
#'  
#' @description This function creates, from a data frame, a list of
#' \code{deadata} objects.
#'
#' @usage make_malmquist(datadea,
#'                nper = NULL,
#'                percol = NULL,
#'                arrangement  = c("horizontal", "vertical"),
#'                ...) 
#'              
#' @param datadea Dataframe with DEA data.
#' @param nper Number of time periods (with dataset in wide format).
#' @param percol Column of time period (with dataset in long format).
#' @param arrangement Horizontal with data in wide format. Vertical with data in long format.
#' @param ... Other options to be passed to the make_deadata function
#'
#' @return An object of class \code{deadata}
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
#' # Example 1. If you have a dataset in wide format.
#' data("Economy")
#' data_example <- make_malmquist(datadea = Economy, 
#'                                nper = 5, 
#'                                arrangement = "horizontal",
#'                                ni = 2, 
#'                                no = 1)
#' # This is the same as:
#' data_example <- make_malmquist(datadea = Economy,
#'                                nper = 5, 
#'                                arrangement = "horizontal",
#'                                inputs = 2:3, 
#'                                outputs = 4)
#' # Example 2. If you have a dataset in long format.
#' data("EconomyLong")
#' data_example2 <- make_malmquist(EconomyLong,
#'                                 percol = 2, 
#'                                 arrangement = "vertical",
#'                                 inputs = 3:4, 
#'                                 outputs = 5)
#'
#' @export

make_malmquist <- function(datadea,
                           nper = NULL,
                           percol = NULL,
                           arrangement  = c("horizontal","vertical"),
                           ...) 
  {
  
  # Checking datadea
  if (!is.data.frame(datadea)) {
    stop("Invalid data datadea (should be a data frame)!")
  }
  datadea <- as.data.frame(datadea)
  
  arrangement <- tolower(arrangement)
  arrangement <- match.arg(arrangement)
  
  if (arrangement == "horizontal") {
    if(is.null(nper)){
      stop("nper must be specified with horizontal arrangement!")
    }
    
    Nio <- ncol(datadea) - 1 # Number of total i/o columns
    NioPeriod <- Nio / nper # Number of i/o columns per period
    datalist <- list()
    for (i in 1:nper){
      datalist[[i]] <- datadea[, c(1, (2 + (i - 1) * NioPeriod):(i * NioPeriod + 1))]
      names(datalist)[i] <- paste("Period",i,sep = ".")
    }
    datalist %>% lapply( function(x) make_deadata(x, ...)) -> dealist
  } else {
    if (is.null(percol)) {
      stop("percol must be supplied with vertical arrangement!")
    }
    periods <- unique(datadea[, percol])
    nper <- length(periods)
    dealist <- list()
    
    for (i in 1:nper) {
      data_temp <- datadea[datadea[, percol] == periods[i], ]
      dealist[[i]] <- make_deadata(datadea = data_temp, ...)
      names(dealist)[i] <- periods[i]
    }
  }
 
  return(dealist)
}