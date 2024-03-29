% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Tone2003}
\alias{Tone2003}
\title{Data: Tone (2003).}
\format{
Data frame with 9 rows and 4 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x}{Input}
  \item{yg}{Output1 ("good" output)}
  \item{yb}{Output2 (undesirable "bad" output)}
}
}
\source{
Tone, K. (2003). "Dealing with undesirable outputs in DEA: A Slacks-Based
Measure (SBM) approach", GRIPS Research Report Series I-2003-0005.
}
\usage{
data("Tone2003")
}
\description{
Data of 9 DMUs producing 2 outputs, being second output undesirable, by using 1 input.
}
\examples{
# Example. Replication of results in Tone (2003), pp 10-11.
data("Tone2003")
data_example <- make_deadata(Tone2003,
                             ni = 1,
                             no = 2,
                             ud_outputs = 2)
result <- model_sbmeff(data_example,
                       rts = "vrs")
efficiencies(result)
targets(result)

}
\seealso{
\code{\link{make_deadata}}, \code{\link{model_sbmeff}}
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
