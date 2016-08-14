make_n_subset <- function(stan_d, n_to_use) {
  # subsets data for model fitting
  indices <- 1:n_to_use
  stan_d$n <- n_to_use
  stan_d$y <- y[indices]
  stan_d$D_site_star <- stan_d$D_site_star[indices, ]
  stan_d$D <- stan_d$D[indices, indices]
  stan_d$X <- as.matrix(stan_d$X[indices, ])
  stan_d
}

get_lp_summary <- function(stanfit, param) {
  summary(stanfit)$summary['lp__', param]
}
