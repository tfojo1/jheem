% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmission.R
\name{set.new.infection.proportions}
\alias{set.new.infection.proportions}
\title{Sets the proportions into which new infections from strata of age x race x subpopulation x sex x risk are distributed}
\usage{
set.new.infection.proportions(jheem, proportions, time = -Inf)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{proportions}{An contact array, with some subset of dimensions: [age.from, race.from, subpopulation.from, sex.from, risk.from, age.to, race.to, subpopulation.to, sex.to, risk.to, non.hiv.subset.to]. from->to follows the direction of HIV transmission}

\item{time}{The time at which these susceptibilities apply, if time-varying}
}
\description{
Sets the proportions into which new infections from strata of age x race x subpopulation x sex x risk are distributed
}
\details{
A new infection proportions array should have dimensions [age, race, subpopulation, sex, risk, non.hiv.subset, continuum, cd4, hiv.subset]. For each combination of age x race x subpopulation x sex x risk x non.hiv.subset, proportions[age,race,subpopulation,sex,risk,non.hiv.subset,,,] should sum to 1 and represent the proportions according to which new infections from the age x race x subpopulation x sex x risk x non.hiv.subset are distributed across strata of continuum x cd4 x hiv.subset
}
\seealso{
Other functions handling new infection proportions: 
\code{\link{get.new.infection.proportions.skeleton}()},
\code{\link{get.uniform.new.infection.proportions}()}

Other functions to specify transmission: 
\code{\link{set.global.transmission.rate}()},
\code{\link{set.susceptibility}()},
\code{\link{set.transmissibility}()},
\code{\link{set.transmission.contact.array}()}
}
\concept{functions handling new infection proportions}
\concept{functions to specify transmission}
