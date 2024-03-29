% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Golany_Roll_1989}
\alias{Golany_Roll_1989}
\title{Data: Golany and Roll (1989).}
\format{
Data frame with 13 rows and 6 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x1}{Input 1}
  \item{x2}{Input 2}
  \item{x3}{Input 3}
  \item{y1}{Output 1}
  \item{y1}{Output 2}
}
}
\source{
Golany, B.; Roll, Y. (1989). "An Application Procedure for DEA". Omega,
International Journal of Management Science, 17(3), 237-250.
\doi{10.1016/0305-0483(89)90029-7}
}
\usage{
data("Golany_Roll_1989")
}
\description{
Data of 13 DMUs using 3 inputs to produce 2 outputs.
}
\examples{
# Example.
data("Golany_Roll_1989")
data_example <- make_deadata(datadea = Golany_Roll_1989,
                             dmus = 1,
                             inputs = 2:4,
                             outputs = 5:6)
result <- cross_efficiency(data_example,
                           orientation = "io",
                           selfapp = TRUE)
result$Arbitrary$cross_eff
result$Arbitrary$e

}
\seealso{
\code{\link{make_deadata}}, \code{\link{model_multiplier}}, \code{\link{cross_efficiency}}
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
