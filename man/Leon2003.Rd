% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Leon2003}
\alias{Leon2003}
\title{Data: Leon, Liern, Ruiz and Sirvent (2003).}
\format{
Data frame with 8 rows and 5 columns. Definition of fuzzy inputs (X)
and fuzzy outputs (Y):
\describe{
  \item{x1}{Input 1}
  \item{alpha1}{spread vector Input 1}
  \item{y1}{Output 1}
  \item{beta1}{spread vector Output 1}
}
}
\source{
Leon, T.; Liern, V. Ruiz, J.; Sirvent, I. (2003). "A Possibilistic
Programming Approach to the Assessment of Efficiency with DEA Models", Fuzzy
Sets and Systems, 139, 407–419. \doi{10.1016/S0165-0114(02)00608-5}
}
\usage{
data("Leon2003")
}
\description{
Data of 8 DMUs with one symmetric triangular fuzzy inputs: Xj = (xj, alphaj), and
one symmetric triangular fuzzy outputs: Yj = (yj, betaj).
}
\examples{
# Example. Replication of results in Leon et. al (2003, p. 416)
data("Leon2003")
data_example <- make_deadata_fuzzy(Leon2003,
                                   dmus = 1,
                                   inputs.mL = 2,
                                   inputs.dL = 3,
                                   outputs.mL = 4,
                                   outputs.dL = 5)
result <- modelfuzzy_possibilistic(data_example,
                                   h = seq(0, 1, by = 0.1),
                                   orientation = "io",
                                   rts = "vrs")
efficiencies(result)

}
\seealso{
\code{\link{make_deadata_fuzzy}}, \code{\link{modelfuzzy_possibilistic}},
\code{\link{cross_efficiency_fuzzy}}, \code{\link{modelfuzzy_guotanaka}}
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
