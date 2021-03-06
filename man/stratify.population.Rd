% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup.R
\name{stratify.population}
\alias{stratify.population}
\title{Substratify a population}
\usage{
stratify.population(population, stratification)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{value}{The value to fill the array with}
}
\description{
Returns an array such that each cell in 'population' is multiplied by the array given by 'stratification' to yield an array of dim=c(dim(population), dim(stratification))
}
\seealso{
\code{\link{get.population.skeleton}}, \code{\link{get.hiv.negative.population.skeleton}}, \code{\link{get.hiv.positive.population.skeleton}}
}
