% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oTools.R
\name{aggRel}
\alias{aggRel}
\title{Calculate aggregate equal weighted reliability}
\usage{
aggRel(data, rel)
}
\arguments{
\item{data}{A dataframe of data}

\item{rel}{An array of reliabilities}
}
\value{
The aggregate reliability
}
\description{
Given a dataframe of data where each column is a measure and each row is a
participant alongside an array of reliabilities (one for each measures),
calculate the aggregate reliability with equal weighting for each measure.
}
