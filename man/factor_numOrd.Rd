% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{factor_numOrd}
\alias{factor_numOrd}
\title{Order factor's levels numerically}
\usage{
factor_numOrd(fct, decreasing = FALSE)
}
\arguments{
\item{fct}{factor. Factor to re-leveled}

\item{decreasing}{logical. Order decreasingly, default = FALSE}
}
\description{
Extracts the first numeric portion of factors' levels to order them numerically
}
\examples{
factor_numOrd(c("100deg", "25deg", "50deg"))
}
