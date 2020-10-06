
#'@export
test <- function()
{
    browser()
}

#'@export
test.dx <- function()
{
    delta = dx_function(hiv_positive = hiv.positive,
                        hiv_negative = hiv.negative,
                        DIMS_HIV = dim(get.hiv.positive.population.skeleton(jheem)),
                        DIMS_NONHIV = dim(get.hiv.negative.population.skeleton(jheem)),
                        AGING_RATES_HIV = pp$AGING_RATES_HIV,
                        AGING_RATES_NONHIV = pp$AGING_RATES_NONHIV,
                        FERTILITY_RATES_HIV = pp$FERTILITY_RATES_HIV,
                        FERTILITY_RATES_NONHIV = pp$FERTILITY_RATES_NONHIV,
                        FRACTION_BIRTHS_INFECTED = pp$FRACTION_BIRTHS_INFECTED,
                        BIRTHS_DIMMASK_NONHIV=parameters$BIRTHS.DIMMASK.NONHIV,
                        BIRTH_PROPORTIONS_NONHIV = pp$BIRTH_PROPORTIONS_NONHIV,
                        BIRTHS_DIMMASK_HIV_TO_NONHIV = parameters$BIRTHS.DIMMASK.HIV.TO.NONHIV,
                        BIRTH_PROPORTIONS_HIV_TO_NONHIV = pp$BIRTH_PROPORTIONS_HIV_TO_NONHIV,
                        BIRTHS_DIMMASK_HIV_TO_HIV = parameters$BIRTHS.DIMMASK.HIV.TO.HIV,
                        BIRTH_PROPORTIONS_HIV_TO_HIV = pp$BIRTH_PROPORTIONS_HIV_TO_HIV,
                        NUM_TRANSMISSION_ROUTES = parameters$NUM_TRANSMISSION_ROUTES,
                        CONTACT_MATRICES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('CONTACT_', i)})],
                        TRANSMISSIBILITIES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('TRANSMISSIBILITY_', i)})],
                        SUSCEPTIBILITIES = pp[sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){paste0('SUSCEPTIBILITY_', i)})],
                        GLOBAL_TRANSMISSION_RATES = sapply(1:parameters$NUM_TRANSMISSION_ROUTES, function(i){pp[[paste0('GLOBAL_TRANSMISSION_RATE_', i)]]}),
                        CONTACT_DIMMASKS = parameters$TRANSMISSION_ROUTE_DIMMASKS,
                        NEW_INFECTION_PROPORTIONS = pp$NEW_INFECTION_PROPORTIONS,
                        HIV_TRANSITION_MATRIX = pp$HIV_POSITIVE_TRANSITION_MATRICES,
                        NONHIV_TRANSITION_MATRIX = pp$HIV_NEGATIVE_TRANSITION_MATRICES,
                        TRACKED_HIV_TRANSITION_INDICES = parameters$TRACKED_HIV_TRANSITION_INDICES-1,
                        TRACKED_HIV_TRANSITION_FROM = parameters$TRACKED_HIV_TRANSITIONS_FROM-1,
                        TRACKED_NONHIV_TRANSITION_INDICES = parameters$TRACKED_NONHIV_TRANSITION_INDICES-1,
                        TRACKED_NONHIV_TRANSITION_FROM = parameters$TRACKED_NONHIV_TRANSITIONS_FROM-1,
                        HIV_SPECIFIC_MORTALITY_RATES = pp$HIV_SPECIFIC_MORTALITY,
                        GENERAL_MORTALITY_RATES_HIV = pp$GENERAL_MORTALITY_FOR_HIV_POSITIVE,
                        GENERAL_MORTALITY_RATES_NONHIV = pp$GENERAL_MORTALITY_FOR_HIV_NEGATIVE,
                        FIXED_SIZE=if (pp$KEEP_STRATA_SIZE_CONSTANT) FIXED_SIZE_BY_STRATA else 0,
                        FIXED_SIZE_STRATA_DIMMASK = parameters$FIXED_SIZE_STRATA_DIMMASK,
                        MODEL_BIRTHS = pp$MODEL_BIRTHS,
                        MODEL_MATERNAL_TRANSMISSION = parameters$MODEL_MATERNAL_TRANSMISSION,
                        TRACK_HIV_SPECIFIC_MORTALITY = parameters$TRACK_HIV_SPECIFIC_MORTALITY,
                        TRACK_HIV_OVERALL_MORTALITY = parameters$TRACK_OVERALL_HIV_MORTALITY,
                        TRACK_NONHIV_MORTALITY = parameters$TRACK_OVERALL_NONHIV_MORTALITY)
}
