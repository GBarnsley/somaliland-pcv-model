function age_values_to_ages(ages)
    cumsum(ages) - ages
end

function age_values_to_age_groups(ages)
    ages_lower = trunc.(Int, age_values_to_ages(ages))
    ages_upper = trunc.(Int, cumsum(ages))
    output = string.(ages_lower) .* "-" .* string.(ages_upper .- 1)
    output[end] = string(ages_lower[end]) * "+"
    output
end

function exponential_proportions(lambda, age_values, unbound_end)
    ages = age_values_to_ages(age_values)

    modelled_props = exp.(-lambda .* ages) .- exp.(-lambda .* (ages .+ age_values))
    if unbound_end
        modelled_props[end] = exp(-lambda * ages[end])
    end
    modelled_props
end

function get_population(age_structure, online)
    #limited info
    #taken from central statistics handbook 2021 see somaliland_population.pdf
    total_population = 4.3*10^6
    #proportions
    #taken from The Somaliland Demographic and Health Survery Table 2.1, see SLHDS2020-Report_2020_table2.1.pdf
    #split into ordered 5 year age groups, proportion of the average household made up by that age group
    population = [5 17.1
    5 16.4
    5 14.9
    5 12.2
    5 7.6
    5 6.2
    5 4.9
    5 4.0
    5 3.6
    5 2.1
    5 3.6
    5 1.7
    5 1.9
    5 0.8
    5 1.3
    5 0.3
    0 1.3]
    #correct final data to make up the last of the age groups in age_structure
    population[end, 1] = sum(age_structure) - sum(population[:, 1])
    #does not sum to 1, likely to due to rounding so correct this
    population[:, 2] .= population[:, 2] ./ sum(population[:, 2])
    #a model requirement is that the population is decreasing
    #hence we'll fit an exponential curve to these values
    #we'll minimize the error with respect to the proportions
    function error_from_props(lambda, population)
        sum(
            (population[1:16, 2] .+ 
                exp.(lambda .* (-1) .* (cumsum(population[1:16, 1]) .- population[1:16, 1])) .* (
                    exp.(lambda .* (-1) .* population[1:16, 1]) .- 1
                )).^2
        ) + ((
            population[17, 2] - exp(-lambda * (sum(population[:, 1]) .- population[17, 1]))
        )^2)
    end
    results = optimize(x -> error_from_props(x[1], population), 1/60, 1)
    exp_rate = Optim.minimizer(results)
    #plot to demonstrate fit to data
    age_group = age_values_to_age_groups(population[:, 1])

    modelled_props = exponential_proportions(exp_rate, population[:, 1], true)

    plot(scatter(1:17, population[:, 2], markerstrokewidth=0, label = "Demographic & Health Survery 2021",
        xticks = (1:17, age_group)))
    plot!(1:17, modelled_props, label = "Fitted Exponential Distribution λ=" * string(trunc(exp_rate, digits = 4)))
    xlabel!("Age Group")
    ylabel!("Proportion of population")
    plot!(size=(900,500))
    savefig("plots/exponential_fit_to_population.png")
    #convert to desired age structure
    new_proportions = exponential_proportions(exp_rate, age_structure, true)
    #convert to population
    final_population = new_proportions .* total_population
    (final_population, exp_rate)
end

function determine_bounds(age1, age1_values, age2, index)
    lower = age1[index]
    upper = age1[index] + age1_values[index]
    lower_age_index = sum(lower .>= age2)
    upper_age_index = sum(upper .>= age2)
    if lower_age_index == upper_age_index
        (true, lower_age_index)
    else
        (false, lower_age_index:(upper_age_index - 1))
    end
end

function get_mixing_matrix(age_structure, online, somaliland_population_rate, somaliland_total_pop)
    if online
        download(
            "https://github.com/kieshaprem/synthetic-contact-matrices/raw/master/generate_synthetic_matrices/output/syntheticmatrices/synthetic_contacts_2021.csv",
            "data/raw/contact_matrix_eth.csv"
        )
        matrices_df = CSV.read(
            "data/raw/contact_matrix_eth.csv",
            DataFrame
        )
        #convert this into a matrix for the given age groups
        filter!(:iso3c => ==("ETH"), matrices_df)
        filter!(:setting => ==("overall"), matrices_df)
        filter!(:location_contact => ==("all"), matrices_df)
        select!(matrices_df, :age_contactor, :age_cotactee, :mean_number_of_contacts)
        matrix = unstack(matrices_df, :age_contactor, :mean_number_of_contacts)
        age_group = string.(matrix[:, 1])
        select!(matrix, age_group)
        writedlm("data/raw/contact_matrix_eth.csv", Array(matrix), ",")
        writedlm("data/raw/contact_matrix_eth_age_groups.csv", age_group, ",")
    end

    #read in age groups and extract lower bound on ages
    ages = readdlm("data/raw/contact_matrix_eth_age_groups.csv")[:, 1]
    ages[end] = parse(Int, first.(ages[end], 2))
    age_values = ages[2:size(ages)[1]] .- ages[1:(size(ages)[1] - 1)]
    push!(age_values, (sum(age_structure) - ages[end]))
    #contactee is on the vertical, values are mean number of contacts a day
    raw_mixing_matrix = readdlm("data/raw/contact_matrix_eth.csv", ',')

    #get ethiopias population data from wpp (https://population.un.org/wpp/)
    eth_pop_raw = readdlm("data/raw/wpp_2022_Ethiopia_ages.csv", ',')
    mixing_age_index = map(x -> (sum(x .>= ages)), eth_pop_raw[:, 1])
    eth_pop = map(x -> sum(eth_pop_raw[mixing_age_index .== x, 2]), 1:(size(ages)[1])) .* 1000

    #assume this captures some level of mixing that is independant on the size of the contactor
    mixing_per_age_per_age = raw_mixing_matrix #TEMP Not sure this holds

    #calculate the somaliland population from the exponential distribution
    sml_pop = exponential_proportions(somaliland_population_rate, age_values, true) .* somaliland_total_pop

    #assume this is constant over the sub age groups and = a weighted average when combining age groups (based on somaliland data)
    target_ages = age_values_to_ages(age_structure)
    n_age = size(age_structure)[1]
    new_mixing = zeros((n_age, n_age))
    for target_index_i in 1:n_age
        index_i_is_sub, index_i_bounds = determine_bounds(target_ages, age_structure, ages, target_index_i)
        for target_index_j in 1:n_age
            index_j_is_sub, index_j_bounds = determine_bounds(target_ages, age_structure, ages, target_index_j)
            if index_i_is_sub & index_j_is_sub
                #if both are sub-age groups
                new_mixing[target_index_i, target_index_j] = mixing_per_age_per_age[index_i_bounds, index_j_bounds]
            elseif !index_i_is_sub & !index_j_is_sub
                #both are larger than one age group
                #take a population weighted average of all values
                new_mixing[target_index_i, target_index_j] = 
                    sum(mixing_per_age_per_age[index_i_bounds, index_j_bounds] .* (sml_pop[index_i_bounds] .* reshape(sml_pop[index_j_bounds], (1, size(index_j_bounds)[1])))) / 
                        sum((sml_pop[index_i_bounds] .* reshape(sml_pop[index_j_bounds], (1, size(index_j_bounds)[1]))))
            elseif index_i_is_sub
                new_mixing[target_index_i, target_index_j] = sum(mixing_per_age_per_age[index_i_bounds, index_j_bounds] .* 
                    sml_pop[index_j_bounds]) / sum(sml_pop[index_j_bounds])
            elseif index_j_is_sub
                new_mixing[target_index_i, target_index_j] = sum(mixing_per_age_per_age[index_i_bounds, index_j_bounds] .* 
                    sml_pop[index_i_bounds]) / sum(sml_pop[index_i_bounds])
            end
        end
    end
    
    #count contacts as going both ways (would be alright if there's not differential on mis classifying ages)
    new_mixing .+= transpose(new_mixing) 

    #rescale so that total level of activity remains the same as the Ethiopia matrix
    scaled_mixing_matrix = new_mixing .* (sum(raw_mixing_matrix) / sum(new_mixing))
    
    #should be fine?
    scaled_mixing_matrix
end