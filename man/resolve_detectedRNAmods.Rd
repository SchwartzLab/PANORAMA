% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{resolve_detectedRNAmods}
\alias{resolve_detectedRNAmods}
\title{Resolve RNAmod multicalling}
\usage{
resolve_detectedRNAmods(insTable, models)
}
\arguments{
\item{insTable}{data.table}

\item{models}{list}
}
\value{
data.table
}
\description{
Using a percentiles version of the linear scores of the conflicting RNAmods
the RNAmod with the highest score according to their percentile will be chosen
and written in the column 'detected_RNAmod_resolve'.
}
\details{
If still more than one RNAmod that reaches the maximum score then those that
reach the maximum score will be output.
}
