library(tidyverse)
library(JuliaCall)
library(drjacoby)

walk(file.path("R", list.files("R")), source)

#get julia functions setup
julia_setup()
julia_source(file.path("julia", "functions.jl"))
julia_source(file.path("julia", "fixed_demography.jl"))

## Prepare Repo Structure
#should we download data from online resources
online <- TRUE

if (!dir.exists(file.path("data", "derived"))){
    dir.create(file.path("data", "derived"))
}

age_structure <- c(
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    c(1, 1, 1, 1, 1, 4, 5, 5, 10, 10, 10, 70) * 12
) #in months

setup_julia_env(online, age_structure)

#define parameters
params <- define_params(
    name = "crate_nvt", min = 0, max = Inf,
    name = "crate_vt", min = 0, max = Inf,
    name = "competition", min = 0, max = 1,
    name = "beta_vt_u5", min = 0, max = Inf,
    name = "beta_vt_u10", min = 0, max = Inf,
    name = "beta_vt_o10", min = 0, max = Inf,
    name = "beta_nvt_u5", min = 0, max = Inf,
    name = "beta_nvt_u10", min = 0, max = Inf,
    name = "beta_nvt_o10", min = 0, max = Inf
)

#define priors + extra params
misc <- list(
    crate_nvt = list(shape = 3, rate = 1/2),
    crate_vt = list(shape = 3, rate = 1/2),
    competition = list(shape1 = 1, shape2 = 1),
    beta_vt_u5 = list(shape = 9, rate = 1/2),
    beta_vt_u10 = list(shape = 9, rate = 1/2),
    beta_vt_o10 = list(shape = 9, rate = 1/2),
    beta_nvt_u5 = list(shape = 9, rate = 1/2),
    beta_nvt_u10 = list(shape = 9, rate = 1/2),
    beta_nvt_o10 = list(shape = 9, rate = 1/2),
    mixing_matrix = julia_eval("mixing_matrix"),
    reps_per_beta = c(
        sum(cumsum(age_structure) < 5*12),
        sum(cumsum(age_structure) < 10*12 & cumsum(age_structure) >= 5*12),
        sum(cumsum(age_structure) >= 10*12)
    ),
    min_prob = 10^-100
)

##load data
#temp for now
data <- c(
    nvt = 800,
    vt = 1000,
    b = 500
)
data <- c(s = julia_eval("sum(N)") - sum(data), data)
data <- list(
    prevalence = data
)
#run drjacoby, can't be in parallel

mcmc_results <- run_mcmc(
    data, params, misc, loglikelihood, logprior,
    burnin = 1000, samples = 1000, rungs = 5, chains = 5
)
