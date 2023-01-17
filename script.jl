using OrdinaryDiffEq, Plots, Turing, Optim

## Load Functions and Model
include("./model/fixed_demography.jl")
include("./functions/functions.jl")


## Prepare Repo Structure

#should we download data from online resources
online = true

if !isdir("data/derived")
    mkdir("data/derived")
end

#sizes of each consecutive age group
age_structure = [1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1/12, 1, 1, 1, 1, 1, 4, 5, 5, 10, 10, 10, 70]

## Load fixed parameters

population = get_population(age_structure, online)
mixing_matrix = get_mixing(age_structure, online)

## Setup Inputs


## Load Data


## Setup Priors


## Fit Model


## Simulate non catch-up roll-out


## Simulate roll-out with catch-up