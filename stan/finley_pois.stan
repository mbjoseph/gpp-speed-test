data {
  int<lower = 1> n;
  int<lower = 1> p;
  matrix[n, p] X;
  int<lower = 1, upper = n> m;
  int<lower = 0> y[n];
  matrix[m, m] D_star;
  matrix[n, m] D_site_star;
}

parameters {
  vector[p] beta;
  real<lower = 0> eta;
  real<lower = 0> sigma;
  real<lower = 0> phi;
  vector[m] w_z;
  vector[n] e_z;
}

transformed parameters {
  vector[n] w;
  vector[n] sigma_e_tilde;
  matrix[m, m] Cstar;
  vector[m] w_star;
  matrix[m, m] inv_Cstar;
  matrix[n, m] C_site_star;
  matrix[n, m] C_ss_inv_Cstar;
  real eta_sq;
  real sig_sq;

  eta_sq = pow(eta, 2);
  sig_sq = pow(sigma, 2);
  
  // latent gp at knots
  for (i in 1:(m-1)) {
    for (j in (i + 1):m) {
      Cstar[i, j] = eta_sq * exp(-D_star[i, j] * phi);
      Cstar[j, i] = Cstar[i, j];
    }
  }

  for (k in 1:m) Cstar[k, k] = eta_sq + sig_sq;
  inv_Cstar = inverse(Cstar);
  w_star = cholesky_decompose(Cstar) * w_z;

  // latent gp at sample locations
  C_site_star = eta_sq * exp(-D_site_star * phi);
  C_ss_inv_Cstar = C_site_star * inv_Cstar;
  w = C_site_star * inv_Cstar * w_star;
  
  // bias adjustment from Finley et al. 2009
  sigma_e_tilde = eta_sq + sig_sq - rows_dot_product(C_ss_inv_Cstar, C_site_star);
  for (i in 1:n) {
    w[i] = w[i] + e_z[i] * sqrt(sigma_e_tilde[i]);
  }
  

}

model {
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  eta ~ normal(0, 1);
  phi ~ normal(0, 5);
  w_z ~ normal(0, 1);
  e_z ~ normal(0, 1);
  y ~ poisson_log(X * beta + w);
}
