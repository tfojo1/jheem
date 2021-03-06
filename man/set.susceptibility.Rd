% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmission.R
\name{set.susceptibility}
\alias{set.susceptibility}
\title{Sets susceptibility to infection for the HIV-negative population}
\usage{
set.susceptibility(
  jheem,
  susceptibility,
  transmission.route.names = jheem$transmission.routes,
  time = -Inf
)
}
\arguments{
\item{jheem}{An initialized JHEEM object}

\item{susceptibility}{An array of susceptibility rates, whose dimensions match the HIV-negative population [age, race, subpopulation, sex, risk, non.hiv.subset].}

\item{transmission.route.names}{The names of the transmission routes to which the transmissibility should apply (a subset of the 'transmission.route.names' given to \code{\link{initialize.jheem}})}

\item{time}{The time at which these susceptibilities apply, if time-varying}
}
\description{
Susceptibility is a characteristic of the HIV-negative population, a multiplier indicating how likely an HIV-negative individual in each stratum is to aquire infection (per partnership) relative to a baseline. ie a value of 1 would not change the susceptibility from the global (baseline) rate
}
\details{
The transmission rate (for one route of transmission) from stratum a (hiv-positive) to stratum b (hiv-negative) is a product of five components (four of which are set as parameters, and one of which - HIV prevalence - is determined by running the model):
\itemize{
 \item {\strong{A global transmission rate} - applies to all strata equally (\code{\link{set.global.transmission.rate}})}
 \item {\strong{Transmissibility} - all factors which depend only on the characteristics of the HIV-positive stratum a (\code{\link{set.transmissibility}})}
 \item {\strong{Contact rate} - all factors which depend on both the general characteristics of the HIV-positive stratum a (ie age, race, subpopulation, sex, risk factor) and the characteristics of the HIV-negative stratum b (\code{\link{set.transmission.contact.array}})}
 \item {\strong{Susceptibility} - all factors which depend only on the characteristics of the HIV-negative stratum b (\code{\link{set.susceptibility}})}
 \item {\strong{Prevalence of HIV} within the general characteristics of the HIV-positive stratum a (ie age, race, subpopulation, sex, risk factor). This is not set, but is estiated at model run time}
}
}
\seealso{
\code{\link{create.susceptibility.array.from.marginals}}, \code{\link{get.hiv.negative.population.skeleton}}

Other functions to specify transmission: 
\code{\link{set.global.transmission.rate}()},
\code{\link{set.new.infection.proportions}()},
\code{\link{set.transmissibility}()},
\code{\link{set.transmission.contact.array}()}
}
\concept{functions to specify transmission}
