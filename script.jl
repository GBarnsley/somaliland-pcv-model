using Zygote, SciMLSensitivity, OrdinaryDiffEq, Plots, Turing, Optim, CSV, DataFrames, DelimitedFiles, LSODA, RecursiveArrayTools

#todo: set up this script, readme, and urban/rural population

## Switch to an AD backend that's compatible with ODE solving
setadbackend(:zygote)

## Load Functions and Model
include("./model/fixed_demography.jl")
include("./functions/functions.jl")

## Prepare Repo Structure

#should we download data from online resources
online = false

if !isdir("data/derived")
    mkdir("data/derived")
end

#sizes of each consecutive age group (in terms of months)
age_structure = cat([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 4, 5, 5, 10, 10, 10, 70] .* 12, dims = 1)

## Load fixed parameters
n_comp = 4
n_age = size(age_structure)[1]
n_dose = 1
n_vacc = 1 + (n_dose*2)

N, somaliland_population_rate = get_population(age_structure, online)
mixing_matrix = get_mixing_matrix(age_structure, online, somaliland_population_rate, sum(N))

serotypes_in_pcv = ()

## Setup Inputs

N = reshape(N, (1, n_age))
age_rate = reshape(1 ./ (age_structure[begin:(end - 1)] .* (365/12)), (1, n_age - 1))

t_simulate = 365*5

t_fit = 365*100
dim_tuple_fitting = (n_comp, n_age, 1)
initial_state = zeros(dim_tuple_fitting)
initial_state[1, :, 1] .= 1
initial_state[1:4, 5, 1] .= 0.25

#temporary vaccination values for fitting
n_dose_fitting = 0
vaccine_efficacy = 0.0
vaccine_waning = 0.0
vaccine_coverage = Float64[0]
vaccine_dose_compartments = Int[0]

## Load Data

#just generate some random for now
target_prevalence = [0.6, 0.2, 0.1, 0.1]

## Setup Priors

#prior_crate_nvt = setup_prior_crate(serotypes_in_pcv, prevalence, true)
#prior_crate_vt = setup_prior_crate(serotypes_in_pcv, prevalence, false)

prior_crate_nvt = 30
prior_crate_vt = 30

prior_competition = [1, 1]

prior_beta_vt_u5 = [9, 0.5]
prior_beta_vt_u10 = [9, 0.5]
prior_beta_vt_o10 = [9, 0.5]
betas_to_age_ns = [17, 2, 5]

prior_beta_nvt_u5 = [9, 0.5]
prior_beta_nvt_u10 = [9, 0.5]
prior_beta_nvt_o10 = [9, 0.5]

## Fit Model
preassigned_arrays = setup_preassigned_arrays(dim_tuple_fitting)

priors = VectorOfArray([
    prior_competition,
    prior_crate_nvt,
    prior_crate_vt,
    prior_beta_vt_u5,
    prior_beta_vt_u10,
    prior_beta_vt_o10,
    prior_beta_nvt_u5,
    prior_beta_nvt_u10,
    prior_beta_nvt_o10,
    betas_to_age_ns,
    mixing_matrix
])

#Cannot use tuples due to AD
fixed_parameters = VectorOfArray([
    collect(dim_tuple_fitting),
    N, 
    age_rate, 
    mixing_matrix, 1/5, mixing_matrix, 1/5, 0.5, #parameters derived from priors (just placeholder values of the same type)
    n_dose_fitting, 
    vaccine_efficacy, 
    vaccine_waning, 
    vaccine_coverage, 
    vaccine_dose_compartments, 
    collect(preassigned_arrays)
])
isconcretetype(VectorOfArray{Any, 2})

pcvm_prob = ODEProblem(pcvm!, initial_state, t_fit, fixed_parameters)

@model function turing_model(target_prevalence, priors, fixed_parameters, pcvm_prob)
    prior_competition, prior_crate_nvt, prior_crate_vt, prior_beta_vt_u5, prior_beta_vt_u10, prior_beta_vt_o10, prior_beta_nvt_u5, prior_beta_nvt_u10, prior_beta_nvt_o10, betas_to_age_ns, mixing_matrix = priors

    #draw from priors
    competition ~ Beta(prior_competition[1], prior_competition[2])

    cdur_nvt ~ Poisson(prior_crate_nvt)
    crate_nvt = 1/cdur_nvt
    cdur_vt ~ Poisson(prior_crate_vt)
    crate_vt = 1/cdur_vt

    beta_vt_u5 ~ Gamma(prior_beta_vt_u5[1], prior_beta_vt_u5[2])
    beta_vt_u10 ~ Gamma(prior_beta_vt_u10[1], prior_beta_vt_u10[2])
    beta_vt_o10 ~ Gamma(prior_beta_vt_o10[1], prior_beta_vt_o10[2])

    beta_vt = cat(repeat([beta_vt_u5], betas_to_age_ns[1]), repeat([beta_vt_u10], betas_to_age_ns[2]), repeat([beta_vt_o10], betas_to_age_ns[3]), dims = 1)' .* 
        mixing_matrix
    
    beta_nvt_u5 ~ Gamma(prior_beta_nvt_u5[1], prior_beta_nvt_u5[2])
    beta_nvt_u10 ~ Gamma(prior_beta_nvt_u10[1], prior_beta_nvt_u10[2])
    beta_nvt_o10 ~ Gamma(prior_beta_nvt_o10[1], prior_beta_nvt_o10[2])

    beta_nvt = cat(repeat([beta_nvt_u5], betas_to_age_ns[1]), repeat([beta_nvt_u10], betas_to_age_ns[2]), repeat([beta_nvt_o10], betas_to_age_ns[3]), dims = 1)' .* 
        mixing_matrix

    #beta_vt_u5 = 0.5
    #beta_vt_u10 = 1
    #beta_vt_o10 = 0.5
    #beta_vt = cat(repeat([beta_vt_u5], betas_to_age_ns[1]), repeat([beta_vt_u10], betas_to_age_ns[2]), repeat([beta_vt_o10], betas_to_age_ns[3]), dims = 1)' .* 
    #    mixing_matrix
    #beta_nvt_u5 = 0.5
    #beta_nvt_u10 = 1
    #beta_nvt_o10 = 0.5
    #beta_nvt = cat(repeat([beta_nvt_u5], betas_to_age_ns[1]), repeat([beta_nvt_u10], betas_to_age_ns[2]), repeat([beta_nvt_o10], betas_to_age_ns[3]), dims = 1)' .* 
    #    mixing_matrix
    #crate_vt = 1/20
    #crate_nvt = 1/20
    #competition = 0.5

    #setup model parameters
    fixed_parameters[4] = beta_nvt
    fixed_parameters[5] = crate_nvt
    fixed_parameters[6] = beta_vt
    fixed_parameters[7] = crate_vt
    fixed_parameters[8] = competition

    #solve ODE
    sol = solve(pcvm_prob, Tsit5(), p = fixed_parameters, maxiters = 10^10)

    #compare to data
    modelled_prev = sum(sol[:, :, :, end] .* N, dims = (2, 3)) ./ sum(N)

    alpha_0 = 1
    alpha = alpha_0 .* modelled_prev[:, 1, 1]

    target_prevalence ~ Dirichlet(alpha)
    
end

sample(turing_model(target_prevalence, priors, fixed_parameters, pcvm_prob), NUTS(), 10)

## Set up vaccination campaigns


## Simulate Vaccination campaigns
