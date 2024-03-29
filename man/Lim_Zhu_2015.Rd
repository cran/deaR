% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Lim_Zhu_2015}
\alias{Lim_Zhu_2015}
\title{Data: Lim and Zhu (2015).}
\format{
Data frame with 37 rows and 7 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x1 = Budget}{Budget}
  \item{y1 = Indirect_economic}{Indirect economic contribution}
  \item{y2 = Direct_economic}{Direct economic contribution}
  \item{y3 = Technical}{Technical contribution}
  \item{y4 = Social}{Social contribution}
  \item{y5 = Scientific}{Scientific contribution}
}
}
\source{
Lim, S.; Zhu, J. (2015). "DEA Cross-Efficiency Under Variable Returns
to Scale". Journal of Operational Research Society, 66(3), p. 476-487.
\doi{10.1057/jors.2014.13}
}
\usage{
data("Lim_Zhu_2015")
}
\description{
Data of 37 R&D project proposal relating to the Turkish iron and steel industry.
Authors consider one input and five outputs.
}
\examples{
# Example. Arbitrary formulation.
# Input-oriented model under variable returns-to-scale.
data("Lim_Zhu_2015")
data_example <- make_deadata(Lim_Zhu_2015,
                             dmus = 1,
                             ni = 1,
                             no = 5)
cross <- cross_efficiency(data_example,
                          epsilon = 0,
                          orientation = "io",
                          rts = "vrs",
                          selfapp = TRUE,
                          M2 = FALSE,
                          M3 = FALSE)
cross$Arbitrary$e

}
\seealso{
\code{\link{make_deadata}}, \code{\link{model_multiplier}},
\code{\link{cross_efficiency}}
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
