logprior <- function(params, misc) {
    dgamma(params["crate_nvt"], misc$crate_nvt$shape, misc$crate_nvt$rate, log = TRUE) +
        dgamma(params["crate_vt"], misc$crate_vt$shape, misc$crate_vt$rate, log = TRUE) +
        dbeta(params["competition"], misc$competition$shape1, misc$competition$shape2, log = TRUE) +
        dgamma(params["beta_vt_u5"], misc$beta_vt_u5$shape, misc$beta_vt_u5$rate, log = TRUE) +
        dgamma(params["beta_vt_u10"], misc$beta_vt_u10$shape, misc$beta_vt_u10$rate, log = TRUE) +
        dgamma(params["beta_vt_o10"], misc$beta_vt_o10$shape, misc$beta_vt_o10$rate, log = TRUE) +
        dgamma(params["beta_nvt_u5"], misc$beta_nvt_u5$shape, misc$beta_nvt_u5$rate, log = TRUE) +
        dgamma(params["beta_nvt_u10"], misc$beta_nvt_u10$shape, misc$beta_nvt_u10$rate, log = TRUE) +
        dgamma(params["beta_nvt_o10"], misc$beta_nvt_o10$shape, misc$beta_nvt_o10$rate, log = TRUE)
}
