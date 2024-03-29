% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transitions.R
\name{set.tracked.transition.hiv.positive}
\alias{set.tracked.transition.hiv.positive}
\title{Specify an HIV-positive transition to track}
\usage{
set.tracked.transition.hiv.positive(
  jheem,
  name,
  transition.dimension,
  from = NULL,
  to = NULL,
  keep.dimensions = c("age", "race", "subpopulation", "sex", "risk", "continuum",
    "cd4", "hiv.subset"),
  age = NULL,
  race = NULL,
  subpopulation = NULL,
  sex = NULL,
  risk = NULL,
  continuum = NULL,
  cd4 = NULL,
  hiv.subset = NULL
)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{name}{The name by which to refer to this tracked transition}

\item{transition.dimension}{The dimension to which this transition applies}

\item{from}{The names of the to and from states within the transition dimension that define the transition}

\item{to}{The names of the to and from states within the transition dimension that define the transition}

\item{keep.dimensions}{Which dimensions to keep in the tracked array (non-keep dimensions are marginalized out)}

\item{age}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{race}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{subpopulation}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{sex}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{risk}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{continuum}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{cd4}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}

\item{hiv.subset}{The subset of the given dimensions for which to track transitions. If null, tracks the transitions across all values of the dimension}
}
\description{
Specify an HIV-positive transition to track
}
