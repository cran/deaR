% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Fortune500}
\alias{Fortune500}
\title{Data: Zhu (2014).}
\format{
Data frame with 15 rows and 6 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x1 = Assets}{Assets (millions of dollars)}
  \item{x2 = Equity}{Equity (millions of dollars)}
  \item{x3 = Employees}{Number of employees}
  \item{y1 = Revenue}{Revenue (millions of dollars)}
  \item{y2 = Profit}{Profit (millions of dollars)}
}
}
\source{
Zhu, J. (2014). Quantitative Models for Performance Evaluation and
Benchmarking. Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer,
New York. \doi{10.1007/978-3-319-06647-9}
}
\usage{
data("Fortune500")
}
\description{
This dataset consists of 15 firms from the Fortune 500 list 1995
(\url{https://fortune.com/ranking/fortune500/}) with 3 inputs and 2 outputs.
}
\examples{
data("Fortune500")
data_Fortune <- make_deadata(datadea = Fortune500,
                             dmus = 1,
                             inputs = 2:4,
                             outputs = 5:6)
result <- model_multiplier(data_Fortune,
                           epsilon = 1e-6,
                           orientation = "io",
                           rts = "crs")
# results for General Motors and Ford Motor are not shown
# by deaR because the solution is infeasible
efficiencies(result)
multipliers(result)

}
\seealso{
\code{\link{make_deadata}}, \code{\link{model_multiplier}}
}
\author{
\strong{Vicente Coll-Serrano} (\email{vicente.coll@uv.es}).
\emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}

\strong{Vicente Bolos} (\email{vicente.bolos@uv.es}).
\emph{Department of Business Mathematics}

\strong{Rafael Benitez} (\email{rafael.suarez@uv.es}).
\emph{Department of Business Mathematics}

University of Valencia (Spain)
}
\keyword{datasets}
