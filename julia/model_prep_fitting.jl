using DifferentialEquations, Optim, CSV, DataFrames, DelimitedFiles, Plots

## Load fixed parameters
n_comp = 4
n_age = size(age_structure)[1]
n_dose = 1
n_vacc = 1 + (n_dose*2)

N, somaliland_population_rate = get_population(age_structure, online);
mixing_matrix = get_mixing_matrix(age_structure, online, somaliland_population_rate, sum(N));

## Setup Inputs

N = reshape(N, (1, n_age));
age_rate = reshape(1 ./ (age_structure[begin:(end - 1)] .* (365/12)), (1, n_age - 1));

t_simulate = 365*5;

t_fit = 365*100;
dim_tuple_fitting = (n_comp, n_age, 1);
initial_state = zeros(dim_tuple_fitting);
initial_state[1, :, 1] .= 1;
initial_state[1:4, 5, 1] .= 0.25;

#temporary vaccination values for fitting
n_dose_fitting = 0;
vaccine_efficacy = 0.0;
vaccine_waning = 0.0;
vaccine_coverage = Float64[0];
vaccine_dose_compartments = Int[0];

## Fit Model
preassigned_arrays = setup_preassigned_arrays(dim_tuple_fitting);

#Cannot use tuples due to AD
fixed_parameters = [
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
];

pcvm_prob = ODEProblem(pcvm!, initial_state, t_fit, fixed_parameters);