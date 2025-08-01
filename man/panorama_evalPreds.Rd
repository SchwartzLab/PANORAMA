% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{panorama_evalPreds}
\alias{panorama_evalPreds}
\title{Evaluate predictions}
\usage{
panorama_evalPreds(
  PANORAMA,
  nucBias = TRUE,
  roundDig = 4,
  removeOtherRNAmods = TRUE,
  rmNotModGenes = TRUE,
  blacklist = NULL,
  summarizeSETS = TRUE
)
}
\arguments{
\item{PANORAMA}{list. PANORAMA object}

\item{nucBias}{logical. Nucleotide bias: only cases in the base nucleotide
of the respective RNAmod are considered.}

\item{roundDig}{numeric. Number of digits to round the resulting specificity,
sensitivity, FDR, and query rate values.}

\item{removeOtherRNAmods}{logical. Remove other RNAmods from training step}
}
\description{
Evaluate predictions using known modified nucleotides identity
}
