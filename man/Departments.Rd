% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sets.R
\docType{data}
\name{Departments}
\alias{Departments}
\title{Data: Tomkins and Green (1988).}
\format{
Data frame with 20 rows and 11 columns. Definition of inputs (X) and outputs (Y):
\describe{
  \item{x1 = Staff}{Average Full Time Academic Staff 82/3-84/5)}
  \item{x2 = Salaries}{1984-5 Salaries Academics and Related (in pounds))}
  \item{x3 = Other_Exp}{1984-5 Other Expenses (in pounds)}
  \item{y1 = Undergrad}{Average Number Undergraduates 82/3-84/5}
  \item{y2 = Research_post}{Research Postgraduates }
  \item{y3 = Taught_post}{Taught Postgraduates}
  \item{y4 = Res_co_income}{Research council income (in pounds)}
  \item{y5 = Other_res_income}{Other research income (in pounds)}
  \item{y6 = Other_income}{Other income (in pounds)}
  \item{y7 = Publications}{Number of publications}
}
}
\source{
Tomkins, C.; Green, R. (1988). "An Experiment in the Use of Data Envelopment
Analysis for Evaluating the Efficiency of UK University Departments of Accounting",
Financial Accountability and Management, 4(2), 147-164.
\doi{10.1111/j.1468-0408.1988.tb00296.x}
}
\usage{
data("Departments")
}
\description{
Data from 20 University accounting departments in the UK.
}
\examples{
# Example.
# Replication of results DEA1 in Tomkins and Green (1988)
data("Departments")
# Calculate Total income
Departments$Total_income <- Departments[, 5] + Departments[, 6] + Departments[, 7]
data_example <- make_deadata(Departments,
                             inputs = 9,
                             outputs = c(2, 3, 4, 12))
result <- model_basic(data_example,
                      orientation = "io",
                      rts = "crs")
efficiencies(result) # Table 3 (p.156)
references(result) # Table 3 (p.157)

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
