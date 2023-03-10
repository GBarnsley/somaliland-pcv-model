loglikelihood <- function(params, data, misc) {
    #generate beta matrices
    beta_vt <- sweep(
        misc$mixing_matrix,
        MARGIN = 2,
        c(
            rep(params["beta_vt_u5"], misc$reps_per_beta[1]),
            rep(params["beta_vt_u10"], misc$reps_per_beta[2]),
            rep(params["beta_vt_o10"], misc$reps_per_beta[3])
        ),
        "*"
    )
    beta_nvt <- sweep(
        misc$mixing_matrix,
        MARGIN = 2,
        c(
            rep(params["beta_nvt_u5"], misc$reps_per_beta[1]),
            rep(params["beta_nvt_u10"], misc$reps_per_beta[2]),
            rep(params["beta_nvt_o10"], misc$reps_per_beta[3])
        ),
        "*"
    )
    #set model parameters in julia
    julia_assign("beta_nvt", beta_nvt)
    julia_assign("crate_nvt", params["crate_nvt"])
    julia_assign("beta_vt", beta_vt)
    julia_assign("crate_vt", params["crate_vt"])
    julia_assign("competition", params["competition"])
    walk(c(
        "fixed_parameters[4] = beta_nvt;",
        "fixed_parameters[5] = crate_nvt;",
        "fixed_parameters[6] = beta_vt;",
        "fixed_parameters[7] = crate_vt;",
        "fixed_parameters[8] = competition;"
    ), julia_command)

    #solve ode in julia
    julia_command(
        "sol = solve(pcvm_prob, p = fixed_parameters);"
    )

    #extract final prevalence
    positive_cases <- julia_eval(
        "sum(sol[:, :, :, end] .* N, dims = (2, 3))[:, 1, 1] ./ sum(N)"
    )
    #ensure greater than 0
    positive_cases[positive_cases < misc$min_prob] <- misc$min_prob

    #multinomial likelihood
    dmultinom(
        x = data$prevalence, prob = positive_cases,
        log = TRUE
    )
}
