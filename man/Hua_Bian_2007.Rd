% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Hua_Bian_2007}
\alias{Hua_Bian_2007}
\title{Data: Hua and Bian (2007).}
\format{
Data frame with 30 rows and 6 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x1 = D-Input1}{Desirable Input 1}
  \item{x2 = D-Input2}{Desirable Input 2}
  \item{y1 = D-Output1}{Desirable Output 1}
  \item{y2 = D-Output2}{Desirable Output 2}
  \item{y3 = UD-Output1}{Undesirable Output 1}
}
}
\source{
Hua Z.; Bian Y. (2007). DEA with Undesirable Factors. In: Zhu J.,
Cook W.D. (eds) Modeling Data Irregularities and Structural Complexities
in Data Envelopment Analysis. Springer, Boston, MA.
\doi{10.1007/978-0-387-71607-7_6}
}
\usage{
data("Hua_Bian_2007")
}
\description{
Data of 30 DMUs with two desirable inputs, two desirable outputs and one udesirable output.
}
\examples{
# Example. Replication of results in Hua and Bian (2007).
data("Hua_Bian_2007")
# The third output is an undesirable output
data_example <- make_deadata(Hua_Bian_2007,
                             ni = 2,
                             no = 3,
                             ud_outputs = 3)

# Translation parameter (vtrans_o) is set to 1500
result <- model_basic(data_example,
                      orientation = "oo",
                      rts = "vrs",
                      vtrans_o = 1500)
eff <- efficiencies(result)
1 / eff # results M5 in Table 6-5 (p.119)

}
\seealso{
\code{\link{make_deadata}}, \code{\link{model_basic}}
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
