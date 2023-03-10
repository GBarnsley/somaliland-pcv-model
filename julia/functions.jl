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
    #lambda, age_values, unbound_end = (exp_rate, population[:, 1], true)
    
    age_values_adjust = age_values ./ 12
    ages = age_values_to_ages(age_values_adjust)

    modelled_props = exp.(-lambda .* ages) .- exp.(-lambda .* (ages .+ age_values_adjust))
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
    population[end, 1] = (sum(age_structure)/12) - sum(population[:, 1])
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

    population[:, 1] .*= 12

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
    #age1, age1_values, age2, index = (target_age_pop, target_age_structure, pop_age, i)
    lower = age1[index]
    upper = age1[index] + age1_values[index]
    lower_age_index = sum(lower .>= age2)
    upper_age_index = sum(upper .>= age2)
    if lower_age_index == upper_age_index || (lower_age_index + 1) == (upper_age_index)
        (true, lower_age_index)
    else
        (false, lower_age_index:(upper_age_index - 1))
    end
end

function restructure_population(raw_pop, pop_age, target_age_structure)
    #raw_pop, pop_age, target_age_structure = (eth_pop_raw[:, 2], eth_pop_raw[:, 1], age_values)
    #raw_pop, pop_age, target_age_structure = (eth_pop_raw[:, 2], eth_pop_raw[:, 1], age_structure)
    
    pop_age_structure = zeros(size(pop_age))
    pop_age_structure[1:(size(pop_age)[1] - 1)] .= pop_age[2:(size(pop_age)[1])] .- pop_age[1:(size(pop_age)[1] - 1)]
    pop_age_structure[end] = sum(target_age_structure) - sum(pop_age_structure)

    target_age_pop = age_values_to_ages(target_age_structure)
    
    new_pop = zeros(size(target_age_structure)[1])
    for i in 1:(size(target_age_structure)[1])
        is_sub, index = determine_bounds(target_age_pop, target_age_structure, pop_age, i)
        if is_sub
            new_pop[i] = raw_pop[index] * target_age_structure[i] / pop_age_structure[index] 
        else
            new_pop[i] = sum(raw_pop[index])
        end
    end
    new_pop
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

    target_ages = age_values_to_ages(age_structure)

    #read in age groups and extract lower bound on ages
    ages = readdlm("data/raw/contact_matrix_eth_age_groups.csv")[:, 1]
    ages[end] = parse(Int, first.(ages[end], 2))
    ages .*= 12 #make months
    age_values = ages[2:size(ages)[1]] .- ages[1:(size(ages)[1] - 1)]
    push!(age_values, (sum(age_structure) - ages[end]))
    #contactee is on the vertical, values are mean number of contacts a day
    raw_mixing_matrix = readdlm("data/raw/contact_matrix_eth.csv", ',')

    #get ethiopias population data from wpp (https://population.un.org/wpp/)
    eth_pop_raw = readdlm("data/raw/wpp_2022_Ethiopia_ages.csv", ',')

    eth_pop = restructure_population(eth_pop_raw[:, 2], eth_pop_raw[:, 1], age_values) .* 1000
    eth_pop_correct_age_structure = restructure_population(eth_pop_raw[:, 2], eth_pop_raw[:, 1], age_structure) .* 1000

    #for computation simplicty we'll normalise these values
    mx_value = maximum([maximum(eth_pop), maximum(eth_pop_correct_age_structure)])
    eth_pop ./= mx_value
    eth_pop_correct_age_structure ./= mx_value
    
    #now convert to the age structure needed
    converted_matrix_partial = zeros((size(age_structure)[1], size(age_values)[1]))
    for i in 1:size(age_structure)[1]
        is_sub, index = determine_bounds(target_ages, age_structure, ages, i)
        for j in 1:size(age_values)[1]
            if is_sub
                #if its just an age group within a larger one we just set the value to be the same
                converted_matrix_partial[i, j] = raw_mixing_matrix[index, j]
            else
                #if there are multiple age groups within it we set the value to their weighted means
                converted_matrix_partial[i, j] = sum(raw_mixing_matrix[index, j] .* eth_pop[index]) ./ sum(eth_pop[index])
            end
        end
    end

    #now expand out in terms of the contactees
    converted_matrix = zeros(size(age_structure)[1], size(age_structure)[1])
    for i in 1:size(age_structure)[1]
        is_sub, index = determine_bounds(target_ages, age_structure, ages, i)
        if is_sub
            #find all indexes of the new age group in this index
            is_sub_other, indexes = determine_bounds(ages, age_values, target_ages, index)
            #now we ensure that our subcompartments sum to the same pop
            multipler = eth_pop_correct_age_structure[i]/sum(eth_pop_correct_age_structure[indexes])
        end
        for j in 1:size(age_structure)[1]
            if is_sub
                #if its just an age group within a larger one we assume that contact rates are uniform within that age group
                converted_matrix[j, i] = converted_matrix_partial[j, index] * multipler
            else
                #if there are multiple age groups within it we set the value to their sum
                converted_matrix[j, i] = sum(converted_matrix_partial[j, index])
            end
        end
    end

    #these should be nearly identical but they are not!
    #sum(converted_matrix_partial .* eth_pop_correct_age_structure, dims = 2)
    #sum(converted_matrix .* eth_pop_correct_age_structure, dims = 2)
    #(sum(raw_mixing_matrix .* eth_pop), sum(converted_matrix_partial .* eth_pop_correct_age_structure), sum(converted_matrix .* eth_pop_correct_age_structure))

    #get the somaliland population
    #sml_pop = exponential_proportions(somaliland_population_rate, age_structure, true) .* somaliland_total_pop

    #convert to somaliland context by assuming that contacts are proportion to age group size as a proportion of the total population
    finalized_matrix = converted_matrix #.* reshape(((sml_pop ./ sum(sml_pop)) ./ (eth_pop_correct_age_structure ./ sum(eth_pop_correct_age_structure))), (1, size(age_structure)[1]))

    finalized_matrix
end