using OrdinaryDiffEq, Plots, Turing, Optim, CSV, DataFrames, DelimitedFiles

#todo: set up this script, readme, and urban/rural population

## Load Functions and Model
include("./model/fixed_demography.jl")
include("./functions/functions.jl")


## Prepare Repo Structure

#should we download data from online resources
online = false

if !isdir("data/derived")
    mkdir("data/derived")
end

#sizes of each consecutive age group
age_structure = [1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1, 1, 1, 1, 1, 4, 5, 5, 10, 10, 10, 70]

## Load fixed parameters
n_comp = 4
n_age = size(age_structure)[1]
n_dose = 1
n_vacc = 1 + (n_dose*2)

N, somaliland_population_rate = get_population(age_structure, online)
mixing_matrix = get_mixing_matrix(age_structure, online, somaliland_population_rate, sum(N))

N = reshape(N, (1, n_age))
age_rate = reshape(1 ./ (age_structure[begin:(end - 1)] .* 365), (1, n_age - 1))

## Setup Inputs


## Load Data


## Setup Priors


## Fit Model


## Simulate non catch-up roll-out


## Simulate roll-out with catch-up