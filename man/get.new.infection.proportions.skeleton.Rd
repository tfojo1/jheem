% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmission.R
\name{get.new.infection.proportions.skeleton}
\alias{get.new.infection.proportions.skeleton}
\title{Get a skeleton new infection proportions array}
\usage{
get.new.infection.proportions.skeleton(jheem, value = 0)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{value}{The value to fill the array with.}
}
\description{
Gets an array, with all values set to value, with dimensions and names corresponding to the new infection proportions array: [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]
}
\details{
A new infection proportions array should have dimensions [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]. For each combination of age x race x subpopulation x sex x risk x non.hiv.subset, proportions[age,race,subpopulation,sex,risk,non.hiv.subset,,,] should sum to 1 and represent the proportions according to which new infections from the age x race x subpopulation x sex x risk x non.hiv.subset are distributed across strata of continuum x cd4 x hiv.subset
}
\seealso{
Other functions handling new infection proportions: 
\code{\link{get.uniform.new.infection.proportions}()},
\code{\link{set.new.infection.proportions}()}
}
\concept{functions handling new infection proportions}
