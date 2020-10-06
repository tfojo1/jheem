
##----------------------------------------------------------------##
##--                         TESTING.R                          --##
##--                                                            --##
##-- A suite of functions to run after internal modifications   --##
##--  to the JHEEM code.                                        --##
##-- These functions compare JHEEM output to output from        --##
##--  simpler systems of ODEs that test specific components     --##
##--  of the model to make sure the output matches.             --##
##-- They are all bundled together in the 'test.jheem' function --##
##--                                                            --##
##----------------------------------------------------------------##



##---------------------------------------------------##
##--        THE PUBLIC-FACING MAIN FUNCTION        --##
##-- Which calls all the individual test functions --##
##---------------------------------------------------##

#'@title Test JHEEM code to make sure that it is working correctly
#'
#'@param error.behavior One of three choices: 'stop' (in which case the test function throws an error and stops), 'browser' (in which case the test function calls the browser), or 'continue' (in which case the function prints an error but continues running)
#'
#'@keywords internal
#'
#'@export
test.jheem <- function(error.behavior=c('stop','browser','continue')[1])
{
    if (!any(c('stop', 'browser', 'continue') == error.behavior))
        stop("error.behavior must be once of 'stop', 'browser', or 'continue'")

    tests = list(test.aging,
                 test.births,
                 test.mortality,
                 test.transitions,
                 test.transmission,
                 test.fixed.population
                 )

    cat("------------------------------------------------------------------\n")
    cat("RUNNING ", length(tests), " TEST FUNCTION", ifelse(length(tests)==1, '', 'S'), "\n", sep='')
    cat("------------------------------------------------------------------\n")

    test.results = sapply(1:length(tests), function(i){
        tests[[i]](i, length(tests), error.behavior)
    })

    cat("------------------------------------------------------------------\n")
    cat(sum(test.results), " OUT OF ", length(test.results), " TESTS RAN WITHOUT ERROR.\n", sep='')
    cat("------------------------------------------------------------------\n")
}

#Tolerance for declaring two results are 'equal'
FRACTIONAL.TOLERANCE = .01

##-----------------------------------##
##-- THE INDIVIDUAL TEST FUNCTIONS --##
##-----------------------------------##

#Aging - everything is fixed except for aging
#uses age strata
test.aging <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Aging Rates...", sep='')

    #Settings
    age.cutoffs = c(0,20,50, Inf)
    pop.per = 100
    start.year = 1
    end.year = 10

    #Set up and run the JHEEM
    jheem = initialize.jheem(age.cutoffs = age.cutoffs, verbose=F, transmission.route.names='sexual')
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per)
    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)

    #Set up the simple ODEs
    age.spans = age.cutoffs[-1] - age.cutoffs[-length(age.cutoffs)]
    n.age = length(age.spans)

    init = rep(pop.per, 2*n.age)
    dx.fn <- function(time, y, params)
    {
        aging.out = y / c(age.spans, age.spans)
        aging.in = c(0, aging.out[1:(n.age-1)],
                     0, aging.out[(n.age+1):(2*n.age-1)])
        list(aging.in - aging.out)
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)

    #compare the results
    result.dim.names = list(1:10, c('hiv.negative', 'hiv.positive'), jheem$age$labels)
    arr.jheem = arr.ode = array(NA, dim=sapply(result.dim.names, length), dimnames=result.dim.names)

    for (age in 1:n.age)
    {
        arr.jheem[,'hiv.negative',age] = access(jheem.results$hiv.negative, age=age)
        arr.jheem[,'hiv.positive',age] = access(jheem.results$hiv.positive, age=age)

        arr.ode[,'hiv.negative',age] = ode.results[1+1:(end.year-start.year+1), 1+age]
        arr.ode[,'hiv.positive',age] = ode.results[1+1:(end.year-start.year+1), 1+n.age+age]
    }

    #Check and respond appropriately
    if (all(fractional.diff(arr.jheem, arr.ode) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations, which change by aging alone, do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
}

#Births - using fertility and birth proportions
#uses age strata and sex strata
#tests continuous time varying
test.births <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Fertility/Births...", sep='')


    #Settings
    age.cutoffs = c(0,20,50, Inf)
    pop.per = 100
    start.year = 1
    end.year = 10

    fertility.rates.by.age.1 = c(0.04, 0.02, 0)
    f.time.1 = -10
    fertility.rates.by.age.2neg = fertility.rates.by.age.1 * 2
    fertility.rates.by.age.2pos = fertility.rates.by.age.1 * 1.8
    f.time.2 = 20

    male.birth.proportion = 0.5

    age.spans = age.cutoffs[-1] - age.cutoffs[-length(age.cutoffs)]
    n.age = length(age.spans)


    #Set up and run the JHEEM
    jheem = initialize.jheem(age.cutoffs = age.cutoffs, sex.strata = c('female','male'), transmission.route.names='sexual', verbose=F)
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per, births = F)

    fertility1 = fertility2pos = fertility2neg = get.general.population.skeleton(jheem, 0)
    for (age in 1:n.age)
    {
        access(fertility1, sex='female', age=age) = fertility.rates.by.age.1[age]
        access(fertility2pos, sex='female', age=age) = fertility.rates.by.age.2pos[age]
        access(fertility2neg, sex='female', age=age) = fertility.rates.by.age.2neg[age]
    }
    jheem = set.fertility(jheem, fertility1, time=f.time.1)
    jheem = set.fertility.hiv.positive(jheem, fertility2pos, time=f.time.2)
    jheem = set.fertility.hiv.negative(jheem, fertility2neg, time=f.time.2)

    birth.prop.arr = get.population.skeleton(jheem, sex=T)
    access(birth.prop.arr, sex='male') = male.birth.proportion
    access(birth.prop.arr, sex='female') = 1 - male.birth.proportion
    birth.proportions = create.birth.proportions.hiv.negative(jheem, birth.prop.arr)
    jheem = set.birth.proportions.no.maternal.transmission(jheem, birth.proportions)

    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)


    #Set up the simple ODEs
    init = rep(pop.per, 2 *n.age * 2)
    dx.fn <- function(time, y, params)
    {
        hiv.pos = array(y[1:(2*n.age)], dim=c(n.age, 2))
        hiv.neg = array(y[2*n.age + 1:(2*n.age)], dim=c(n.age, 2))

        #Aging
        change.pos = -hiv.pos / age.spans
        change.pos[-1,] = change.pos[-1,] - change.pos[-n.age,]

        change.neg = -hiv.neg / age.spans
        change.neg[-1,] = change.neg[-1,] - change.neg[-n.age,]

        #interpolate fertility
        f.pos = fertility.rates.by.age.1 + (fertility.rates.by.age.2pos - fertility.rates.by.age.1) / (f.time.2 - f.time.1) * (time - f.time.1)
        f.neg = fertility.rates.by.age.1 + (fertility.rates.by.age.2neg - fertility.rates.by.age.1) / (f.time.2 - f.time.1) * (time - f.time.1)

        #births
        births = sum(hiv.pos[,1] * f.pos + hiv.neg[,1] * f.neg)
        change.neg[1,2] = change.neg[1,2] + births * male.birth.proportion
        change.neg[1,1] = change.neg[1,1] + births * (1 - male.birth.proportion)

        #return
        list(c(as.numeric(change.pos), as.numeric(change.neg)))
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)


    #Compare the results
    result.dim.names = list(1:10, jheem$age$labels, c('female','male'))
    arr.jheem.pos = arr.ode.pos = arr.jheem.neg = arr.ode.neg =
        array(NA, dim=sapply(result.dim.names, length), dimnames=result.dim.names)

    for (age in 1:n.age)
    {
        arr.jheem.neg[,age,'female'] = access(jheem.results$hiv.negative, age=age, sex='female')
        arr.jheem.neg[,age,'male'] = access(jheem.results$hiv.negative, age=age, sex='male')

        arr.jheem.pos[,age,'female'] = access(jheem.results$hiv.positive, age=age, sex='female')
        arr.jheem.pos[,age,'male'] = access(jheem.results$hiv.positive, age=age, sex='male')
    }

    for (year in 1:10)
    {
        arr.ode.pos[year,,] = ode.results[1+year,1+1:(2*n.age)]
        arr.ode.neg[year,,] = ode.results[1+year,1+2*n.age + 1:(2*n.age)]
    }

    #Check and respond appropriately
    if (all(fractional.diff(arr.jheem.neg, arr.ode.neg) <= FRACTIONAL.TOLERANCE) &&
        all(fractional.diff(arr.jheem.pos, arr.ode.pos) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations, which change by aging and new births, do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
}

#Mortality - both general and HIV-specific
#Uses race, subpopulation, non.hiv.subset, hiv.subset, cd4
#Also tests time varying
test.mortality <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Mortality...", sep='')



    #Settings
    pop.per = 100
    start.year = 1
    end.year = 10

    races = c('black','white','hispanic')
    n.race = length(races)
    subpopulations = c('here','there')
    n.subpop = length(subpopulations)
    non.hiv.subsets = c('healthy','sick')
    n.nhsub = length(non.hiv.subsets)
    hiv.subsets = c('healthyH', 'sickH')
    n.hsub = length(hiv.subsets)
    cd4.strata = c('nonAIDS','AIDS')
    n.cd4 = length(cd4.strata)

    general.mortality.rates.by.race = c(.06, .04, .07)
    subpopulation.mortality.multipliers = c(1, 1.1)
    sick.increment = 0.02
    aids.increment.1 = 0.2
    aids.increment.2 = 0.1
    t1 = 2
    t2 = 5

    #Set up and run the JHEEM
    jheem = initialize.jheem(race.strata=races,
                             subpopulation = subpopulations,
                             nonhiv.subsets = non.hiv.subsets,
                             hiv.subsets = hiv.subsets,
                             cd4.strata = cd4.strata,
                             transmission.route.names='sexual',
                             verbose=F)
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per, mortality = F)

    gen.mort = get.general.population.skeleton(jheem, 0)
    for (race in 1:n.race)
        access(gen.mort, race=race) = general.mortality.rates.by.race[race]

    for (subpop in 1:n.subpop)
        access(gen.mort, subpopulation = subpop) = access(gen.mort, subpopulation = subpop) * subpopulation.mortality.multipliers[subpop]

    gen.mort.nonhiv = expand.population.to.hiv.negative(jheem, gen.mort)
    access(gen.mort.nonhiv, non.hiv.subset = 2) = access(gen.mort.nonhiv, non.hiv.subset = 2) + sick.increment
    jheem = set.general.mortality.hiv.negative(jheem, gen.mort.nonhiv)

    gen.mort.hiv = expand.population.to.hiv.positive(jheem, gen.mort)
    access(gen.mort.hiv, hiv.subset = 2) = access(gen.mort.hiv, hiv.subset = 2) + sick.increment
    jheem = set.general.mortality.hiv.positive(jheem, gen.mort.hiv)

    hiv.mort.1 = hiv.mort.2 = get.hiv.positive.population.skeleton(jheem, 0)
    access(hiv.mort.1, cd4=2) = aids.increment.1
    access(hiv.mort.2, cd4=2) = aids.increment.2
    jheem = set.hiv.specific.mortality(jheem, hiv.mort.1, time=t1)
    jheem = set.hiv.specific.mortality(jheem, hiv.mort.2, time=t2)

    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)


    #Set up the simple ODEs
    num.hiv = n.race * n.subpop * n.cd4 * n.hsub
    num.nonhiv = n.race * n.subpop * n.nhsub
    init = rep(pop.per, num.hiv + num.nonhiv)

    non.hiv.mort = array(0, dim=c(n.race, n.subpop, n.nhsub))
    hiv.mort.1 = array(0, dim=c(n.race, n.subpop, n.cd4, n.hsub))
    for (race in 1:n.race)
    {
        non.hiv.mort[race,,] = general.mortality.rates.by.race[race]
        hiv.mort.1[race,,,] = general.mortality.rates.by.race[race]
    }

    for (subpop in 1:n.subpop)
    {
        non.hiv.mort[,subpop,] = non.hiv.mort[,subpop,] * subpopulation.mortality.multipliers[subpop]
        hiv.mort.1[,subpop,,] = hiv.mort.1[,subpop,,] * subpopulation.mortality.multipliers[subpop]
    }

    non.hiv.mort[,,2] = non.hiv.mort[,,2] + sick.increment
    hiv.mort.1[,,,2] = hiv.mort.1[,,,2] + sick.increment

    hiv.mort.2 = hiv.mort.1
    hiv.mort.1[,,2,] = hiv.mort.1[,,2,] + aids.increment.1
    hiv.mort.2[,,2,] = hiv.mort.2[,,2,] + aids.increment.2

    dx.fn <- function(time, y, params)
    {
        hiv.pos = array(y[1:num.hiv], dim=c(n.race, n.subpop, n.cd4, n.hsub))
        hiv.neg = array(y[num.hiv + 1:num.nonhiv], dim=c(n.race, n.subpop, n.nhsub))

        #HIV-negative mortality
        change.neg = -hiv.neg * non.hiv.mort

        #interpolate hiv mortality
        if (time <= t1)
            hiv.mort = hiv.mort.1
        else if (time >= t2)
            hiv.mort = hiv.mort.2
        else
            hiv.mort = hiv.mort.1 + (hiv.mort.2 - hiv.mort.1) / (t2 - t1) * (time-t1)

        #HIV-positive mortality
        change.pos = -hiv.pos * hiv.mort

        #return
        list(c(as.numeric(change.pos), as.numeric(change.neg)))
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)


    #Compare the results
    hiv.result.dim.names = list(1:10, races, subpopulations, cd4.strata, hiv.subsets)
    nonhiv.result.dim.names = list(1:10, races, subpopulations, non.hiv.subsets)

    arr.jheem.pos = arr.ode.pos = array(NA, dim=sapply(hiv.result.dim.names, length), dimnames=hiv.result.dim.names)
    arr.jheem.neg = arr.ode.neg = array(NA, dim=sapply(nonhiv.result.dim.names, length), dimnames=nonhiv.result.dim.names)


    for (year in 1:10)
    {
        for (race in 1:n.race)
        {
            for (subpop in 1:n.subpop)
            {
                arr.jheem.neg[year, race, subpop, ] = access(jheem.results$hiv.negative, year=year, race=race, subpopulation = subpop)
                for (cd4 in 1:n.cd4)
                    arr.jheem.pos[year, race, subpop, cd4, ] = access(jheem.results$hiv.positive, year=year, race=race, subpopulation = subpop, cd4=cd4)
            }
        }

        arr.ode.pos[year,,,,] = ode.results[1+year, 1+1:num.hiv]
        arr.ode.neg[year,,,] = ode.results[1+year, 1+num.hiv+1:num.nonhiv]
    }

    #Check and respond appropriately
    if (all(fractional.diff(arr.jheem.neg, arr.ode.neg) <= FRACTIONAL.TOLERANCE) &&
        all(fractional.diff(arr.jheem.pos, arr.ode.pos) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations, which change only by mortality (general and HIV-specific), do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
}

#also tests transmission, transit
test.fixed.population <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Fixed Population Strata Sizes...", sep='')

    #-- Settings --#
    pop.per = 100
    start.year = 1
    end.year = 10

    age.cutoffs = c(0,20,50, Inf)
    age.spans = age.cutoffs[-1] - age.cutoffs[-length(age.cutoffs)]
    n.age = length(age.spans)
    races = c('white','black')
    n.race = length(races)
    risks = c('non_IDU', 'IDU')
    n.risk = length(risks)
    continuum = c('undiagnosed', 'diagnosed')
    n.cont = length(continuum)

    WHITE=1
    BLACK=2
    NONIDU=1
    IDU=2
    UNDX = 1
    DX = 2

    t.rate = 0.04
    idu.incidence = 0.005
    dx.rate = 0.2
    mortality.by.age = c(.01, .02, .04)
    fertility.rate = 0.04
    fix.pop.until.year = 5

    #-- Set up and run the JHEEM --#
    jheem = initialize.jheem(age.cutoffs = age.cutoffs,
                             race.strata=races,
                             risk.strata = risks,
                             continuum.of.care.states = continuum,
                             transmission.route.names='sexual',
                             verbose=F)
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per,
                           transmission = F,
                           transitions = F,
                           mortality = F,
                           births = F)

    #Transmission
    jheem = set.global.transmission.rate(jheem, t.rate)
    jheem = set.transmissibility(jheem, 1)
    jheem = set.susceptibility(jheem, 1)
    jheem = set.transmission.contact.array(jheem, 1)
    prop = get.uniform.new.infection.proportions(jheem, initial.continuum=1)
    jheem = set.new.infection.proportions(jheem, prop)

    #Transisions (IDU)
    trans.idu = get.transition.array.skeleton(jheem, risk = T)
    access.transition(trans.idu, risk.from='non_IDU', risk.to='IDU') = idu.incidence
    jheem = set.transition.array.hiv.negative(jheem, trans.idu)

    trans.cont = get.transition.array.skeleton(jheem, continuum = T)
    access.transition(trans.cont, continuum.from='undiagnosed', continuum.to='diagnosed') = dx.rate
    trans.hiv = create.hiv.positive.transition.array.from.marginals(jheem, trans.cont, trans.idu)
    jheem = set.transition.array.hiv.positive(jheem, trans.hiv)

    #Mortality
    mort = get.general.population.skeleton(jheem)
    for (age in 1:n.age)
        access(mort, age=age) = mortality.by.age[age]
    jheem = set.general.mortality(jheem, mort)
    jheem = set.hiv.specific.mortality(jheem, 0)

    #Fertility
    jheem = set.fertility(jheem, fertility.rate)
    prop = create.birth.proportions.hiv.negative(jheem, risk.born.into.state = 1)
    jheem = set.birth.proportions.no.maternal.transmission(jheem, prop)

    #Fix pop
    jheem = set.fixed.population.strata(jheem, fix.age=T, fix.race=T)
    jheem = set.fixed.population.size(jheem)
    jheem = set.use.fixed.population(jheem, T, -Inf)
    jheem = set.use.fixed.population(jheem, F, fix.pop.until.year)

    #run it
    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)


    #-- Set up the simple ODEs --#
    num.hiv = n.age * n.race * n.risk * n.cont
    num.nonhiv = n.age * n.race * n.risk
    init = rep(pop.per, num.hiv + num.nonhiv)
    init.pos = array(init[1:num.hiv], dim=c(n.age, n.race, n.risk, n.cont))
    init.neg = array(init[num.hiv + 1:num.nonhiv], dim=c(n.age, n.race, n.risk))
    target.pop.size = rowSums(init.pos, dims=2) + rowSums(init.neg, dims=2)

    dx.fn <- function(time, y, params)
    {
        hiv.pos = array(y[1:num.hiv], dim=c(n.age, n.race, n.risk, n.cont))
        hiv.neg = array(y[num.hiv + 1:num.nonhiv], dim=c(n.age, n.race, n.risk))
        change.neg = array(0, dim=dim(hiv.neg))
        change.pos = array(0, dim=dim(hiv.pos))

        #Aging
        hiv.aging.up = hiv.pos / age.spans
        nonhiv.aging.up = hiv.neg / age.spans

        change.pos = change.pos - hiv.aging.up
        change.neg = change.neg - nonhiv.aging.up
        change.pos[-1,,,] = change.pos[-1,,,] + hiv.aging.up[-n.age,,,]
        change.neg[-1,,] = change.neg[-1,,] + nonhiv.aging.up[-n.age,,]


        #Mortality
        hiv.mort = hiv.pos * mortality.by.age
        nonhiv.mort = hiv.neg * mortality.by.age

        change.pos = change.pos - hiv.mort
        change.neg = change.neg - nonhiv.mort


        #HIV transmission
        hiv.prev = sum(hiv.pos) / (sum(hiv.neg) + sum(hiv.pos))
        new.cases.from = hiv.neg * t.rate * hiv.prev

        change.neg = change.neg - new.cases.from
        change.pos[,,,1] = change.pos[,,,1] + new.cases.from


        #IDU transitions
        nonidu.to.idu.hiv = hiv.pos[,,NONIDU,] * idu.incidence
        nonidu.to.idu.nonhiv = hiv.neg[,,NONIDU] * idu.incidence

        change.pos[,,NONIDU,] = change.pos[,,NONIDU,] - nonidu.to.idu.hiv
        change.pos[,,IDU,] = change.pos[,,IDU,] + nonidu.to.idu.hiv

        change.neg[,,NONIDU] = change.neg[,,NONIDU] - nonidu.to.idu.nonhiv
        change.neg[,,IDU] = change.neg[,,IDU] + nonidu.to.idu.nonhiv


        #Continuum transitions
        diagnosis = hiv.pos[,,,UNDX] * dx.rate

        change.pos[,,,UNDX] = change.pos[,,,UNDX] - diagnosis
        change.pos[,,,DX] = change.pos[,,,DX] + diagnosis


        #Births
        if (time >= fix.pop.until.year)
        {
            pop.by.race = apply(hiv.pos, 2, sum) + apply(hiv.neg, 2, sum)
            births = pop.by.race * fertility.rate
            change.neg[1,,1] = change.neg[1,,1] + births
        }

        #Fixed Population
        if (time < fix.pop.until.year)
        {
            unconstrained.target.pos = hiv.pos + change.pos
            unconstrained.target.neg = hiv.neg + change.neg
            unconstrained.pop.size = rowSums(unconstrained.target.pos, dims=2) + rowSums(unconstrained.target.neg, dims=2)

            ratio = target.pop.size / unconstrained.pop.size

            constrained.target.pos = array(NA, dim=dim(hiv.pos))
            constrained.target.neg = array(NA, dim=dim(hiv.neg))

            for (age in 1:n.age)
            {
                for (race in 1:n.race)
                {
                    constrained.target.pos[age,race,,] = unconstrained.target.pos[age,race,,] * ratio[age,race]
                    constrained.target.neg[age,race,] = unconstrained.target.neg[age,race,] * ratio[age,race]
                }
            }

            change.pos[-1,,,] = constrained.target.pos[-1,,,] - hiv.pos[-1,,,]
            change.neg[-1,,] = constrained.target.neg[-1,,] - hiv.neg[-1,,]

            changes.to.first.age = rowSums(change.pos[1,,,]) + rowSums(change.neg[1,,])
            change.neg[1,,1] = change.neg[1,,1] - changes.to.first.age
        }

        #return
        list(c(as.numeric(change.pos), as.numeric(change.neg)))
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)


    #Compare the results
    hiv.result.dim.names = list(1:10, jheem$age$labels, races, risks, continuum)
    nonhiv.result.dim.names = list(1:10, jheem$age$labels, races, risks)

    arr.jheem.pos = arr.ode.pos = array(NA, dim=sapply(hiv.result.dim.names, length), dimnames=hiv.result.dim.names)
    arr.jheem.neg = arr.ode.neg = array(NA, dim=sapply(nonhiv.result.dim.names, length), dimnames=nonhiv.result.dim.names)


    for (year in 1:10)
    {
        for (age in 1:n.age)
        {
            for (race in 1:n.race)
            {
                for (risk in 1:n.risk)
                {
                    arr.jheem.neg[year, age, race, risk] = access(jheem.results$hiv.negative, year=year, age=age, race=race, risk=risk)

                    for (cont in 1:n.cont)
                        arr.jheem.pos[year, age, race, risk, cont] = access(jheem.results$hiv.positive, year=year, age=age, race=race, risk=risk, continuum=cont)
                }
            }
        }

        arr.ode.pos[year,,,,] = ode.results[1+year, 1+1:num.hiv]
        arr.ode.neg[year,,,] = ode.results[1+year, 1+num.hiv+1:num.nonhiv]
    }

    jheem.pop.sums.by.year = rowSums(arr.jheem.neg, dims=3) + rowSums(arr.jheem.pos, dims=3)

    #Check and respond appropriately
    if (!all(apply(jheem.pop.sums.by.year[1:(fix.pop.until.year-1),,], 1, function(pop){all(round(pop)==round(target.pop.size))})))
    {
        error.msg = 'The population sizes did not stay fixed'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
    else if (all(fractional.diff(arr.jheem.neg, arr.ode.neg) <= FRACTIONAL.TOLERANCE) &&
        all(fractional.diff(arr.jheem.pos, arr.ode.pos) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
}

test.transmission <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Transmission...", sep='')

    #-- Settings --#
    pop.per = 100
    start.year = 1
    end.year = 10

    races = c('white','black')
    n.race = length(races)
    sexes = c('female', 'male')
    n.sex = length(sexes)
    risks = c('non_IDU', 'IDU')
    n.risk = length(risks)
    continuum = c('undiagnosed', 'unsuppressed', 'suppressed')
    n.cont = length(continuum)
    non.hiv.subsets = c('no_PrEP', 'PrEP')
    n.nhsub = length(non.hiv.subsets)

    WHITE=1
    BLACK=2
    FEMALE=1
    MALE=2
    NONIDU=1
    IDU=2
    NO_PREP=1
    PREP=2
    UNDX = 1
    UNSUPP = 2
    SUPP = 3

    t.rate = 0.005
    diagnosed.transmission.mult = 0.9
    prep.mult = 0.04

    sex.frequency.by.race = c(10,20)
    black.black.prop = 0.9
    white.white.prop = 0.7
    male.male.prop = 0.2

    idu.frequency = 5

    #-- Set up and run the JHEEM --#
    jheem = initialize.jheem(race.strata=races,
                             sex.strata = sexes,
                             nonhiv.subsets = non.hiv.subsets,
                             risk.strata = risks,
                             continuum.of.care.states = continuum,
                             transmission.route.names = c('sexual','IDU'),
                             verbose=F)
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per, transmission = F)

    #Global rate
    jheem = set.global.transmission.rate(jheem, t.rate)

    #Transmissibility
    trans = get.hiv.positive.population.skeleton(jheem, 0)
    access(trans, continuum='undiagnosed') = 1
    access(trans, continuum='unsuppressed') = diagnosed.transmission.mult
    jheem = set.transmissibility(jheem, trans)

    #Susceptibility
    susc = get.hiv.negative.population.skeleton(jheem, 1)
    access(susc, non.hiv.subset = 'PrEP') = prep.mult
    jheem = set.susceptibility(jheem, susc)

    #Sexual Contact
    sex.contact = get.contact.array.skeleton(jheem, sex=T)
    access(sex.contact, sex.from='male', sex.to='female') = 1
    access(sex.contact, sex.from='female', sex.to='male') = 1-male.male.prop
    access(sex.contact, sex.from='male', sex.to='male') = male.male.prop

    race.contact = get.contact.array.skeleton(jheem, race=T)
    access(race.contact, race.from='black', race.to='black') = black.black.prop
    access(race.contact, race.from='white', race.to='black') = 1-black.black.prop
    access(race.contact, race.from='white', race.to='white') = white.white.prop
    access(race.contact, race.from='black', race.to='white') = 1-white.white.prop

    sex.freq = array(sex.frequency.by.race, dim=n.race, dimnames=list(races))

    contact = create.contact.array.from.marginals(jheem, sex.contact, race.contact, sex.freq)
    jheem = set.transmission.contact.array(jheem, transmission.route.names = 'sexual', contact)

    #IDU contact
    idu.contact = get.contact.array.skeleton(jheem, value=0, risk=T)
    access(idu.contact, risk.from='IDU', risk.to='IDU') = idu.frequency
    jheem = set.transmission.contact.array(jheem, transmission.route.names = 'IDU', idu.contact)

    #Infection proportions
    prop = get.uniform.new.infection.proportions(jheem, initial.continuum=1)
    jheem = set.new.infection.proportions(jheem, prop)

    #run it
    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)


    #-- Set up the simple ODEs --#
    num.hiv = n.race * n.sex * n.risk * n.cont
    num.nonhiv = n.race * n.sex * n.risk * n.nhsub
    init = rep(pop.per, num.hiv + num.nonhiv)

    dx.fn <- function(time, y, params)
    {
        hiv.pos = array(y[1:num.hiv], dim=c(n.race, n.sex, n.risk, n.cont))
        hiv.neg = array(y[num.hiv + 1:num.nonhiv], dim=c(n.race, n.sex, n.risk, n.nhsub))
        change.neg = array(0, dim=dim(hiv.neg))
        change.pos = array(0, dim=dim(hiv.pos))

        #prevalences in strata of race x sex x risk
        hiv.pos.marg = rowSums(hiv.pos, dims=3)
        hiv.neg.marg = rowSums(hiv.neg, dims=3)

        hiv.pos.trans = hiv.pos
        hiv.pos.trans[,,,SUPP] = 0
        hiv.pos.trans[,,,UNSUPP] = hiv.pos.trans[,,,UNSUPP] * diagnosed.transmission.mult
        hiv.pos.trans = rowSums(hiv.pos.trans, dims=3)

        pop = hiv.pos.marg + hiv.neg.marg
        hiv.prev = hiv.pos.trans / pop
        hiv.prev.idu = sum(hiv.pos.trans[,,IDU]) / sum(pop[,,IDU])

        #sexual transmission
        white.male.new = t.rate * sex.frequency.by.race[WHITE] *
            (male.male.prop * white.white.prop * hiv.prev[WHITE, MALE, ] +
                 (1-male.male.prop) * white.white.prop * hiv.prev[WHITE, FEMALE, ] +
                 male.male.prop * (1-white.white.prop) * hiv.prev[BLACK, MALE, ] +
                 (1-male.male.prop) * (1-white.white.prop) * hiv.prev[BLACK, FEMALE, ]
            )
        white.female.new = t.rate * sex.frequency.by.race[WHITE] *
            (white.white.prop * hiv.prev[WHITE, MALE, ] +
                 (1-white.white.prop) * hiv.prev[BLACK, MALE, ]
            )
        black.male.new = t.rate * sex.frequency.by.race[BLACK] *
            (male.male.prop * black.black.prop * hiv.prev[BLACK, MALE, ] +
                 (1-male.male.prop) * black.black.prop * hiv.prev[BLACK, FEMALE, ] +
                 male.male.prop * (1-black.black.prop) * hiv.prev[WHITE, MALE, ] +
                 (1-male.male.prop) * (1-black.black.prop) * hiv.prev[WHITE, FEMALE, ]
            )
        black.female.new = t.rate * sex.frequency.by.race[BLACK] *
            (black.black.prop * hiv.prev[BLACK, MALE, ] +
                 (1-black.black.prop) * hiv.prev[WHITE, MALE, ]
            )

        #* These are each 2-vectors of non-idu, idu
        white.male.new.prep = white.male.new * prep.mult * hiv.neg[WHITE, MALE,, PREP]
        white.male.new.noprep = white.male.new * hiv.neg[WHITE, MALE,, NO_PREP]
        black.male.new.prep = black.male.new * prep.mult * hiv.neg[BLACK, MALE,, PREP]
        black.male.new.noprep = black.male.new * hiv.neg[BLACK, MALE,, NO_PREP]
        white.female.new.prep = white.female.new * prep.mult * hiv.neg[WHITE, FEMALE,, PREP]
        white.female.new.noprep = white.female.new * hiv.neg[WHITE, FEMALE,, NO_PREP]
        black.female.new.prep = black.female.new * prep.mult * hiv.neg[BLACK, FEMALE,, PREP]
        black.female.new.noprep = black.female.new * hiv.neg[BLACK, FEMALE,, NO_PREP]

        #IDU transmission
        idu.new.noprep = t.rate * idu.frequency * hiv.prev.idu * hiv.neg[,,IDU,NO_PREP]
        idu.new.prep = t.rate * idu.frequency * prep.mult * hiv.prev.idu * hiv.neg[,,IDU,PREP]

        #outflows
        change.neg[WHITE, MALE,, PREP] = -white.male.new.prep
        change.neg[WHITE, MALE,, NO_PREP] = -white.male.new.noprep
        change.neg[BLACK, MALE,, PREP] = -black.male.new.prep
        change.neg[BLACK, MALE,, NO_PREP] = -black.male.new.noprep
        change.neg[WHITE, FEMALE,, PREP] = -white.female.new.prep
        change.neg[WHITE, FEMALE,, NO_PREP] = -white.female.new.noprep
        change.neg[BLACK, FEMALE,, PREP] = -black.female.new.prep
        change.neg[BLACK, FEMALE,, NO_PREP] = -black.female.new.noprep

        change.neg[,,IDU,NO_PREP] = change.neg[,,IDU,NO_PREP] - idu.new.noprep
        change.neg[,,IDU,PREP] = change.neg[,,IDU,PREP] - idu.new.prep

        #inflows
        change.pos[WHITE, MALE,, 1] = white.male.new.prep + white.male.new.noprep
        change.pos[BLACK, MALE,, 1] = black.male.new.prep + black.male.new.noprep
        change.pos[WHITE, FEMALE,, 1] = white.female.new.prep + white.female.new.noprep
        change.pos[BLACK, FEMALE,, 1] = black.female.new.prep + black.female.new.noprep

        change.pos[,,IDU,1] = change.pos[,,IDU,1] + idu.new.prep + idu.new.noprep

        #return
        list(c(as.numeric(change.pos), as.numeric(change.neg)))
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)


    #Compare the results
    hiv.result.dim.names = list(1:10, races, sexes, risks, continuum)
    nonhiv.result.dim.names = list(1:10, races, sexes, risks, non.hiv.subsets)

    arr.jheem.pos = arr.ode.pos = array(NA, dim=sapply(hiv.result.dim.names, length), dimnames=hiv.result.dim.names)
    arr.jheem.neg = arr.ode.neg = array(NA, dim=sapply(nonhiv.result.dim.names, length), dimnames=nonhiv.result.dim.names)


    for (year in 1:10)
    {
        for (race in 1:n.race)
        {
            for (sex in 1:n.sex)
            {
                for (risk in 1:n.risk)
                {
                    for (sub in 1:n.nhsub)
                        arr.jheem.neg[year, race, sex, risk, sub] = access(jheem.results$hiv.negative, year=year, race=race, sex=sex, risk=risk, non.hiv.subset=sub)

                    for (cont in 1:n.cont)
                        arr.jheem.pos[year, race, sex, risk, cont] = access(jheem.results$hiv.positive, year=year, race=race, sex=sex, risk=risk, continuum=cont)
                }
            }
        }

        arr.ode.pos[year,,,,] = ode.results[1+year, 1+1:num.hiv]
        arr.ode.neg[year,,,,] = ode.results[1+year, 1+num.hiv+1:num.nonhiv]
    }

    #Check and respond appropriately
    if (all(fractional.diff(arr.jheem.neg, arr.ode.neg) <= FRACTIONAL.TOLERANCE) &&
        all(fractional.diff(arr.jheem.pos, arr.ode.pos) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations, which change only by HIV transmission, do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }

}

#Test transitions
#Uses race, risk, non.hiv.subset, continuum of care
test.transitions <- function(test.num=1, num.tests=1, error.behavior='stop')
{
    cat("Test ", test.num, " of ", num.tests, ": Transitions...", sep='')

    #-- Settings --#
    pop.per = 100
    start.year = 1
    end.year = 10

    races = c('black','white','hispanic')
    n.race = length(races)
    risks = c('non_IDU', 'IDU')
    n.risk = length(risks)
    continuum = c('undiagnosed', 'unsuppressed', 'suppressed')
    n.cont = length(continuum)
    non.hiv.subsets = c('low_risk', 'high_risk')
    n.nhsub = length(non.hiv.subsets)

    NONIDU=1
    IDU=2
    LOW_RISK=1
    HIGH_RISK=2
    UNDX = 1
    UNSUPP = 2
    SUPP = 3

    diagnosis.rate.by.race = c(2, 3, 1.5) / 10
    suppression.rate.by.race = c(1, 2, 1) / 10
    unsuppression.rate.by.race = c(.5, .25, .33) / 10

    incident.idu.rate.hiv = .04
    incident.idu.rate.nonhiv.by.subset = c(.01, .08)
    idu.remission.rate = .2

    #-- Set up and run the JHEEM --#
    jheem = initialize.jheem(race.strata=races,
                             nonhiv.subsets = non.hiv.subsets,
                             risk.strata = risks,
                             continuum.of.care.states = continuum,
                             transmission.route.names='sexual',
                             verbose=F)
    jheem = set.null.jheem(jheem, init.pop.size.per = pop.per, transitions = F)

    #set this one explicitly
    non.hiv.trans = get.transition.array.skeleton(jheem, race=T, risk=T, non.hiv.subset = T)
    for (sub in 1:n.nhsub)
        access.transition(non.hiv.trans, non.hiv.subset = sub, risk.from='non_IDU', risk.to='IDU') = incident.idu.rate.nonhiv.by.subset[sub]
    access.transition(non.hiv.trans, risk.from='IDU', risk.to='non_IDU') = idu.remission.rate
    jheem = set.transition.array.hiv.negative(jheem, non.hiv.trans)

    #set this one by multiplying marginals
    hiv.trans.sub.risk = get.transition.array.skeleton(jheem, risk=T)
    access(hiv.trans.sub.risk, risk.from='non_IDU', risk.to='IDU') = incident.idu.rate.hiv
    access(hiv.trans.sub.risk, risk.from='IDU', risk.to='non_IDU') = idu.remission.rate

    hiv.trans.sub.continuum = array(0, dim=c(n.race, n.cont, n.cont), dimnames=list(races, continuum, continuum))
    for (race in 1:n.race)
    {
        hiv.trans.sub.continuum[race, 1, 2] = diagnosis.rate.by.race[race]
        hiv.trans.sub.continuum[race, 2, 3] = suppression.rate.by.race[race]
        hiv.trans.sub.continuum[race, 3, 2] = unsuppression.rate.by.race[race]
    }

    hiv.trans = create.hiv.positive.transition.array.from.marginals(jheem, hiv.trans.sub.risk, hiv.trans.sub.continuum)
    jheem = set.transition.array.hiv.positive(jheem, hiv.trans)

    #run it
    jheem.results = run.jheem(jheem, start.year, end.year, verbose = F)


    #-- Set up the simple ODEs --#
    num.hiv = n.race * n.risk * n.cont
    num.nonhiv = n.race * n.risk * n.nhsub
    init = rep(pop.per, num.hiv + num.nonhiv)

    dx.fn <- function(time, y, params)
    {
        hiv.pos = array(y[1:num.hiv], dim=c(n.race, n.risk, n.cont))
        hiv.neg = array(y[num.hiv + 1:num.nonhiv], dim=c(n.race, n.risk, n.nhsub))
        change.neg = array(0, dim=dim(hiv.neg))
        change.pos = array(0, dim=dim(hiv.pos))

        #HIV-negative - IDU transitions
        nonidu.to.idu.low.risk = hiv.neg[,NONIDU,LOW_RISK] * incident.idu.rate.nonhiv.by.subset[1]
        nonidu.to.idu.high.risk = hiv.neg[,NONIDU,HIGH_RISK] * incident.idu.rate.nonhiv.by.subset[2]
        idu.to.nonidu = hiv.neg[,IDU,] * idu.remission.rate

        change.neg[,NONIDU,LOW_RISK] = change.neg[,NONIDU,LOW_RISK] - nonidu.to.idu.low.risk
        change.neg[,IDU,LOW_RISK] = change.neg[,IDU,LOW_RISK] + nonidu.to.idu.low.risk

        change.neg[,NONIDU,HIGH_RISK] = change.neg[,NONIDU,HIGH_RISK] - nonidu.to.idu.high.risk
        change.neg[,IDU,HIGH_RISK] = change.neg[,IDU,HIGH_RISK] + nonidu.to.idu.high.risk

        change.neg[,IDU,] = change.neg[,IDU,] - idu.to.nonidu
        change.neg[,NONIDU,] = change.neg[,NONIDU,] + idu.to.nonidu


        #HIV-positive IDU
        nonidu.to.idu = hiv.pos[,NONIDU,] * incident.idu.rate.hiv
        idu.to.nonidu = hiv.pos[,IDU,] * idu.remission.rate

        change.pos[,IDU,] = change.pos[,IDU,] + nonidu.to.idu - idu.to.nonidu
        change.pos[,NONIDU,] = change.pos[,NONIDU,] + idu.to.nonidu - nonidu.to.idu


        #HIV-positive continuum
        for (race in 1:n.race)
        {
            dx.to.unsupp = hiv.pos[race,,UNDX] * diagnosis.rate.by.race[race]
            unsupp.to.supp = hiv.pos[race,,UNSUPP] * suppression.rate.by.race[race]
            supp.to.unsupp = hiv.pos[race,,SUPP] * unsuppression.rate.by.race[race]

            change.pos[race,,UNDX] = change.pos[race,,UNDX] - dx.to.unsupp
            change.pos[race,,UNSUPP] = change.pos[race,,UNSUPP] + dx.to.unsupp + supp.to.unsupp - unsupp.to.supp
            change.pos[race,,SUPP] = change.pos[race,,SUPP] + unsupp.to.supp - supp.to.unsupp
        }

        #return
        list(c(as.numeric(change.pos), as.numeric(change.neg)))
    }
    ode.results = deSolve::ode(y=init, times=start.year:(end.year+1), func=dx.fn, parms=NULL)


    #Compare the results
    hiv.result.dim.names = list(1:10, races, risks, continuum)
    nonhiv.result.dim.names = list(1:10, races, risks, non.hiv.subsets)

    arr.jheem.pos = arr.ode.pos = array(NA, dim=sapply(hiv.result.dim.names, length), dimnames=hiv.result.dim.names)
    arr.jheem.neg = arr.ode.neg = array(NA, dim=sapply(nonhiv.result.dim.names, length), dimnames=nonhiv.result.dim.names)


    for (year in 1:10)
    {
        for (race in 1:n.race)
        {
            for (risk in 1:n.risk)
            {
                for (sub in 1:n.nhsub)
                    arr.jheem.neg[year, race, risk, sub] = access(jheem.results$hiv.negative, year=year, race=race, risk=risk, non.hiv.subset=sub)

                for (cont in 1:n.cont)
                    arr.jheem.pos[year, race, risk, cont] = access(jheem.results$hiv.positive, year=year, race=race, risk=risk, continuum=cont)
            }
        }

        arr.ode.pos[year,,,] = ode.results[1+year, 1+1:num.hiv]
        arr.ode.neg[year,,,] = ode.results[1+year, 1+num.hiv+1:num.nonhiv]
    }

    #Check and respond appropriately
    if (all(fractional.diff(arr.jheem.neg, arr.ode.neg) <= FRACTIONAL.TOLERANCE) &&
        all(fractional.diff(arr.jheem.pos, arr.ode.pos) <= FRACTIONAL.TOLERANCE))
    {
        #print success and return
        cat('OK\n')
        T
    }
    else
    {
        error.msg = 'The populations, which change only by transitions for risk (IDU) and continuum of care, do not match in the JHEEM simulation and the simple ODE simulation'
        if (error.behavior == 'stop')
            stop(error.msg)
        else if (error.behavior == 'browser')
        {
            cat('ERROR!:\n  ', error.msg, '\n', sep='')
            browser()
        }
        else
            cat('ERROR!:\n  ', error.msg, '\n', sep='')

        F
    }
}

##-------------##
##-- HELPERS --##
##-------------##

#A helper to set up JHEEM with null values
# for each param, if true, sets the param to a null value
# if false, does not set that param
set.null.jheem <- function(jheem,
                           init.pop.size.per=100,
                           init.population=T,
                           births=T,
                           mortality=T,
                           transmission=T,
                           transitions=T)
{
    #Init population
    if (init.population)
    {
        init.hiv.negative = get.hiv.negative.population.skeleton(jheem, init.pop.size.per)
        init.hiv.positive = get.hiv.positive.population.skeleton(jheem, init.pop.size.per)
        jheem = set.initial.populations(jheem, init.hiv.negative, init.hiv.positive)
    }

    #Fertility/births
    if (births)
    {
        fertility = get.general.population.skeleton(jheem, 0)
        jheem = set.fertility(jheem, fertility)

        birth.prop = create.birth.proportions.hiv.negative(jheem,
                                                           subpopulations.born.into.state = 1,
                                                           sex.born.into.state = 1,
                                                           risk.born.into.state = 1,
                                                           non.hiv.subset.born.into.state = 1)
        jheem = set.birth.proportions.no.maternal.transmission(jheem, birth.prop)
    }

    #Mortality
    if (mortality)
    {
        mortality = get.general.population.skeleton(jheem, 0)
        jheem = set.general.mortality(jheem, mortality)

        hiv.mortality = get.hiv.negative.population.skeleton(jheem, 0)
        jheem = set.hiv.specific.mortality(jheem, mortality)
    }

    #Transmission
    if (transmission)
    {
        jheem = set.global.transmission.rate(jheem, 0)
        jheem = set.transmissibility(jheem, 1)
        jheem = set.susceptibility(jheem, 1)
        jheem = set.transmission.contact.array(jheem, 0)
    }

    #Transitions
    if (transitions)
    {
        jheem = set.transition.array.hiv.negative(jheem, 0)
        jheem = set.transition.array.hiv.positive(jheem, 0)
    }

    #Return it
    jheem
}

fractional.diff <- function(x, y)
{
    abs(x-y)/x
}
