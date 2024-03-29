% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.friends.R
\name{is.friends}
\alias{is.friends}
\title{Friends check.}
\usage{
is.friends(datadea,
             dmu_eval = NULL,
             dmu_ref = NULL,
             rts = c("crs", "vrs", "nirs", "ndrs"),
             tol = 1e-6)
}
\arguments{
\item{datadea}{The data, including \code{n} DMUs, \code{m} inputs and \code{s} outputs.}

\item{dmu_eval}{A numeric vector containing the subset of DMUs to be checked.
If \code{NULL} (default), all DMUs are considered.}

\item{dmu_ref}{A numeric vector containing which DMUs are the evaluation reference set.
If \code{NULL} (default), all DMUs are considered.}

\item{rts}{A string, determining the type of returns to scale, equal to "crs" (constant),
"vrs" (variable), "nirs" (non-increasing) or "ndrs" (non-decreasing).}

\item{tol}{Numeric, a tolerance margin for checking efficiency. It is 1e-6 by default.}
}
\value{
Returns \code{TRUE} if \code{dmu_eval} is friends of \code{dmu_ref},
        and \code{FALSE} otherwise.
}
\description{
Checks whether a subset of DMUs is friends or not, according to Tone (2010).
}
\examples{
data("PFT1981")
datadea <- make_deadata(PFT1981,
                        ni = 5,
                        no = 3)
subset1 <- c(15, 16, 17, 19) # Subset of DMUs to be checked
result1 <- is.friends(datadea = datadea,
                      dmu_eval = subset1,
                      dmu_ref = 1:20) # We only consider a cluster formed by the first 20 DMUs
subset2 <- c(15, 16, 17, 20) # Another subset of DMUs to be checked
result2 <- is.friends(datadea = datadea,
                      dmu_eval = subset2,
                      dmu_ref = 1:20) # We only consider a cluster formed by the first 20 DMUs
 
}
\references{
Tone, K. (2010). "Variations on the theme of slacks-based measure of efficiency
in DEA", European Journal of Operational Research, 200, 901-907.
\doi{10.1016/j.ejor.2009.01.027}
}
\seealso{
\code{\link{maximal_friends}}, \code{\link{model_sbmeff}}
}
\author{
\strong{Vicente Coll-Serrano} (\email{vicente.coll@uv.es}).
\emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}

\strong{Vicente Bolós} (\email{vicente.bolos@uv.es}).
\emph{Department of Business Mathematics}

\strong{Rafael Benítez} (\email{rafael.suarez@uv.es}).
\emph{Department of Business Mathematics}

University of Valencia (Spain)
}
