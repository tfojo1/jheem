% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmission.R
\name{get.uniform.new.infection.proportions}
\alias{get.uniform.new.infection.proportions}
\title{Get a new infection proportions array where new infections are seeded uniformly}
\usage{
get.uniform.new.infection.proportions(
  jheem,
  initial.continuum = jheem$continuum.of.care[1],
  initial.cd4 = jheem$cd4.strata[1],
  initial.hiv.subset = jheem$hiv.subsets[1]
)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{initial.continuum, initial.cd4, initial.hiv.subset}{The index (indices) or name(s) of the dimensions into which new infections are to be distributed uniformly}
}
\description{
The return value will be a new infection proportions array such that the new infection proportions are the same across all strata of age x race x subpopulation x sex x risk x non.hiv.subset, with new infections being distributed uniformly across the specified continuum, cd4, and hiv.subset specified
}
\details{
A new infection proportions array should have dimensions [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]. For each combination of age x race x subpopulation x sex x risk x non.hiv.subset, proportions[age,race,subpopulation,sex,risk,non.hiv.subset,,,] should sum to 1 and represent the proportions according to which new infections from the age x race x subpopulation x sex x risk x non.hiv.subset are distributed across strata of continuum x cd4 x hiv.subset
}
\seealso{
Other functions handling new infection proportions: 
\code{\link{get.new.infection.proportions.skeleton}()},
\code{\link{set.new.infection.proportions}()}
}
\concept{functions handling new infection proportions}
