setup_julia_env <- function(online, age_structure) {
    julia_command(
        paste0(
            "age_structure = cat([",
            paste0(age_structure, collapse = ", "),
            "], dims = 1);"
        )
    )
    if (online){
        julia_command(
            "online = true;"
        )
    } else {
        julia_command(
            "online = false;"
        )
    }
    julia_source(file.path("julia", "model_prep_fitting.jl"))
}
