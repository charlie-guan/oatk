% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{buildX}
\alias{buildX}
\title{buildX}
\usage{
buildX(dat, mut, ps)
}
\arguments{
\item{dat}{Raw data from the HIV Stanfrod Database}

\item{mut}{list of mutation to keep}

\item{ps}{position of the mutations}
}
\value{
Design matrix X
}
\description{
Helper function to create the design matrix X with the input mutations/positions.
Taken from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R
}
