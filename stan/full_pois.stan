data {
  int<lower = 1> n;
  int<lower = 1> p;
  matrix[n, p] X;
  int<lower = 0> y[n];
  matrix[n, n] D;
}

parameters {
  vector[n] z;
  real<lower=0> eta;
  real<lower=0> phi;
  real<lower=0> sigma;
  vector[p] beta;
}

transformed parameters {
  cov_matrix[n] Sigma;
  vector[n] w;
  real<lower = 0> eta_sq;
  real<lower = 0> sig_sq;
  
  eta_sq = pow(eta, 2);
  sig_sq = pow(sigma, 2);
  
  for (i in 1:(n-1)) {
    for (j in (i + 1):n) {
      Sigma[i, j] = eta_sq * exp(-D[i, j] * phi);
      Sigma[j, i] = Sigma[i, j];
    }
  }

  for (k in 1:n) Sigma[k, k] = eta_sq + sig_sq;
  w = cholesky_decompose(Sigma) * z;
}

model {
  eta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  phi ~ normal(0, 5);
  beta ~ normal(0, 1);
  z ~ normal(0, 1);
  y ~ poisson_log(X * beta + w);
}
