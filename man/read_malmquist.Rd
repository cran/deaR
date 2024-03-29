% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_malmquist.R
\name{read_malmquist}
\alias{read_malmquist}
\title{read_malmquist}
\usage{
read_malmquist(datadea,
               nper = NULL,
               percol = NULL,
               arrangement  = c("horizontal", "vertical"),
               ...)
}
\arguments{
\item{datadea}{Data frame with DEA data.}

\item{nper}{Number of time periods (with dataset in wide format).}

\item{percol}{Column of time period (with dataset in long format).}

\item{arrangement}{Horizontal with data in wide format. Vertical with data in long format.}

\item{...}{Other options to be passed to the \code{make_deadata} function.}
}
\description{
This function is deprecated. Use \code{make_malmquist} instead.
}
