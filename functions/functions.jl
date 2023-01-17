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
    ages = trunc.(Int, cumsum(population[:, 1]) - population[:, 1])
    ages_upper = trunc.(Int, cumsum(population[:, 1]))
    age_group = string.(ages) .* "-" .* string.(ages_upper .- 1)
    age_group[17] = "80+"

    modelled_props = exp.(-exp_rate .* ages) .- exp.(-exp_rate .* ages_upper)
    modelled_props[17] = exp(-exp_rate * ages[17])

    plot(scatter(1:17, population[:, 2], markerstrokewidth=0, label = "Demographic & Health Survery 2021",
        xticks = (1:17, age_group)))
    plot!(1:17, modelled_props, label = "Fitted Exponential Distribution λ=" * string(trunc(exp_rate, digits = 4)))
    xlabel!("Age Group")
    ylabel!("Proportion of population")
    plot!(size=(900,500))
    savefig("plots/exponential_fit_to_population.png")
    #convert to desired age structure
    target_ages = cumsum(age_structure)
    new_proportions = exp.(-exp_rate .* (target_ages .- age_structure)) .- exp.(-exp_rate .* target_ages)
    #convert to population
    final_population = (new_proportions ./ sum(new_proportions)) .* total_population
    #final check
    #sum(diff(final_population .* 1 ./ age_structure) .< 0) == (size(age_structure)[1] - 1)
    final_population
end