
function setup_preassigned_arrays(dim_tuple)
    #N, aging, births, migrants, force_of_infection_nvt, force_of_infection_vt, s_to_nvt, vt_to_b, s_to_vt, nvt_to_b, nvt_to_s, vt_to_s, b_to_vt, b_to_nvt
    n_comp, n_age, n_vacc = dim_tuple
    (
        zeros((n_age, 1)), #force_of_infection_nvt
        zeros((n_age, n_vacc)), #force_of_infection_vt
        zeros((n_age, n_vacc)), #s_to_nvt
        zeros((n_age, n_vacc)), #vt_to_b
        zeros((n_age, n_vacc)), #s_to_vt
        zeros((n_age, n_vacc)), #nvt_to_b
        zeros((n_age, n_vacc)), #nvt_to_s
        zeros((n_age, n_vacc)), #vt_to_s
        zeros((n_age, n_vacc)), #b_to_vt
        zeros((n_age, n_vacc)), #b_to_nvt
        1:(n_age - 1), #age_offset_1
        2:n_age, #age_offset_2
        zeros((n_comp, n_age - 1, n_vacc)), #vaccination_u
        zeros((1, n_age - 1, 1)), #vaccinated
        zeros((1, n_age - 1, 1)), #unvaccinated
        zeros((n_comp, n_age, n_vacc)) #waning
    )
end

function calculate_force_of_infection!(foi, inf, beta, vaccine_efficacy, n_age)
    #foi, inf, beta, vaccine_efficacy, n_age = (force_of_infection_vt, VT .+ B, beta_vt, vaccine_efficacy, n_age)
    #this is not efficent and needs optimisation
    foi .= sum(beta .* reshape(sum(inf, dims = 2), (1, n_age)), dims = 2) .* (1 .- vaccine_efficacy)
end

#! as this modifies du in place, saving on memory
function pcvm!(du, u, p, t)
    #get parameters and compartments, we don't name them otherwise as its quicker
    dim_tuple, N, age_rate, beta_nvt, crate_nvt, beta_vt, crate_vt, competition, n_dose, vaccine_efficacy, vaccine_waning, vaccine_coverage, vaccine_dose_compartments, preassigned_arrays = p
    n_comp, n_age, n_vacc = dim_tuple
    force_of_infection_nvt, force_of_infection_vt, s_to_nvt, vt_to_b, s_to_vt, nvt_to_b, nvt_to_s, vt_to_s, b_to_vt, b_to_nvt, age_offset_1, age_offset_2, vaccination_u, vaccinated, unvaccinated, waning = preassigned_arrays
    dS = @view du[1, :, :]
    dNVT = @view du[2, :, :]
    dVT = @view du[3, :, :]
    dB = @view du[4, :, :]
    S = @view u[1, :, :]
    NVT = @view u[2, :, :]
    VT = @view u[3, :, :]
    B = @view u[4, :, :]
    #vaccinations (based upon ageing so this calculates before the ageing)
    vaccination_u .= u[:, age_offset_1, :] #to hold the effective proportions accounting for vaccinations
    if n_dose > 0
        for dose in 1:n_dose

            vaccinated .= sum(u[:, :, vaccine_dose_compartments[dose + 1]], dims = (1, 3))[:, age_offset_1, :] 
            unvaccinated .= sum(u[:, :, vaccine_dose_compartments[dose]], dims = (1, 3))[:, age_offset_1, :]

            #cannot preassing this as it is of variable length
            to_be_vaccinated = (((vaccine_coverage[:, age_offset_1, dose] .* (vaccinated .+ unvaccinated)) .- vaccinated) ./ unvaccinated) .* u[:, age_offset_1, vaccine_dose_compartments[dose]]

            #does not apply when coverage is higher than the target
            to_be_vaccinated[to_be_vaccinated .< 0] .= 0.0
            to_be_vaccinated[isnan.(to_be_vaccinated)] .= 0.0

            vaccination_u[:, :, vaccine_dose_compartments[dose + 1][1]] .+= sum(to_be_vaccinated, dims = 3)
            vaccination_u[:, :, vaccine_dose_compartments[dose]] .-= to_be_vaccinated
        end
    end
    #ageing and deaths + births
    older_age_groups = @view u[:, age_offset_2, :]
    births = age_rate[1, end] .* N[1, n_age - 1] #assume that its the same as deaths in the final compartment
    du[:, 1, :] .= (-births ./ N[:, 1]) .* u[:, 1, :]
    du[1, 1, 1] = (births / N[1, 1]) * (1 - u[1, 1, 1])
    du[:, age_offset_2, :] .= age_rate .* (N[:, age_offset_1] ./ N[:, age_offset_2]) .* (vaccination_u .- older_age_groups)
    #waning (can also model a delay in protection if setup correctly)
    waning .= u .* vaccine_waning
    du .-= waning
    du[:, :, 2:n_vacc] .+= waning[:, :, 1:(n_vacc - 1)]
    #transmission
    #this is frequency dependent since we're calculating with proportions
    ##NVT
    calculate_force_of_infection!(force_of_infection_nvt, NVT .+ B, beta_nvt, 0.0, n_age)
    s_to_nvt .= S .* force_of_infection_nvt
    vt_to_b .= VT .* force_of_infection_nvt .* competition
    ##VT
    calculate_force_of_infection!(force_of_infection_vt, VT .+ B, beta_vt, vaccine_efficacy, n_age)
    s_to_vt .= S .* force_of_infection_vt
    nvt_to_b .= NVT .* force_of_infection_vt .* competition
    ##clearance
    nvt_to_s .= NVT .* crate_nvt
    vt_to_s .= VT .* crate_vt
    b_to_vt .= B .* crate_nvt
    b_to_nvt .= B .* crate_vt
    dS .+= nvt_to_s .+ vt_to_s .- s_to_vt .- s_to_nvt
    dNVT .+= b_to_nvt .+ s_to_nvt .- nvt_to_b .- nvt_to_s
    dVT .+= b_to_vt .+ s_to_vt .- vt_to_b .- vt_to_s
    dB .+= nvt_to_b .+ vt_to_b .- b_to_nvt .- b_to_vt
end
