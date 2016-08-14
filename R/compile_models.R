library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

m_inits <- list()
for (m in models) {
  m_inits[m] <- stan(m, data = stan_d, chains = 1, iter = 1, par = 'lp__')
}
