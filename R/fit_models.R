# Comparing MCMC efficiency: GPP vs. GP across a gradient of sample size
library(dplyr)
library(gridExtra)

# make data, compile models, and load helper functions -------------------------
nmax <- 300 # max. number of sites
nvec <- floor(exp(seq(log(20), log(nmax), length.out = 7)))

models <- list.files('stan', '.stan', full.names = TRUE)

res <- expand.grid(n = nvec, model = models, 
                       elapsed_t = NA, n_eff = NA, Rhat = NA)

source('R/generate_data.R')
source('R/zzz.R')
source('R/compile_models.R')

# set up simulation experiment, containers for results, mcmc params ------------
n_trials <- nrow(res)
fits <- vector(length = n_trials, mode = 'list')
pars_to_watch <- c('eta', 'sigma', 'phi', 'beta', 'w')
n_it <- 1000
n_chains <- 2

# iterate over models and sample sizes, storing estimates, n_eff, and Rhat -----
for (i in 1:n_trials) {
  res$elapsed_t[i] <- system.time(
    fits[[i]] <- stan(fit = m_inits[[res$model[i]]], 
                      iter = n_it, chains = n_chains, 
                      data = make_n_subset(stan_d, res$n[i]), 
                      par = pars_to_watch)
  )['elapsed']
}


# get effective samples sizes, Rhat, & MCMC efficiency -------------------------
res$n_eff <- sapply(fits, get_lp_summary, param = 'n_eff')
res$Rhat <- sapply(fits, get_lp_summary, param = 'Rhat')
res <- res %>%
  mutate(efficiency = n_eff / elapsed_t, 
         Model = ifelse(grepl('finley', model), 
                        'GPP', 'Full GP'))

save.image()
