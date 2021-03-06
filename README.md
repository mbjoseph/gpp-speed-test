# Gaussian predictive process models in Stan
Max Joseph  
August 14, 2016  

Gaussian process (GP) models are computationally demanding for large datasets.
Much work has been done to circumvent expensive matrix operations that arise in parameter estimation with larger datasets via sparse and/or reduced rank covariance matrices ([Datta et al. 2016](http://arxiv.org/pdf/1406.7343.pdf) provide a nice review).
What follows is an implementation of a spatial Gaussian predictive process Poisson GLM in [Stan](http://mc-stan.org/), following [Finley et al. 2009](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2743161/): *Improving the performance of predictive process modeling for large datasets*, with a comparison to a full rank GP model in terms of execution time and MCMC efficiency.

## Generative model

Assume we have $n$ spatially referenced counts $y(s)$ made at spatial locations $s_1, s_2, ..., s_n$, which depend on a latent mean zero Gaussian process with an isotropic stationary exponential covariance function:

$$y(s) \sim \text{Poisson}(\text{exp}(X^T(s) \beta + w(s)))$$

$$w(s) \sim GP(0, C(d))$$

$$[C(d)]_{i, j} = \eta^2 \text{exp}(- d_{ij} \phi)) + I(i = j) \sigma^2$$

where $y(s)$ is the response at location $s$, $X$ is an $n \times p$ design matrix, $\beta$ is a length $p$ parameter vector, $\sigma^2$ is a "nugget" parameter, $\eta^2$ is the variance parameter of the Gaussian process, $d_{ij}$ is a spatial distance between locations $s_i$ and $s_j$, and $\phi$ determines how quickly the correlation in $w$ decays as distance increases.
For simplicity, the point locations are uniformly distributed in a 2d unit square spatial region.



To estimate this model in a Bayesian context, we might be faced with taking a Cholesky decomposition of the $n \times n$ matrix $C(d)$ at every iteration in an MCMC algorithm, a costly operation for large $n$.

## Gaussian predictive process representation

Computational benefits of Gaussian predictive process models arise from the estimation of the latent Gaussian process at $m << n$ locations (knots).
Instead of taking the Cholesky factorization of the $n \times n$ covariance matrix, we instead factorize the $m \times m$ covariance matrix corresponding to the covariance in the latent spatial process among knots. 
Knot placement is a non-trivial topic, but for the purpose of illustration let's place knots on a grid over the spatial region of interest. 
Note that here the knots are fixed, but it is possible to model knot locations stochastically as in [Guhaniyogi et al. 2012](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268014/).

![](README_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

Knots are shown as stars and the points are observations.

We will replace $w$ above with an estimate of $w$ that is derived from a reduced rank representation of the latent spatial process.
Below, the vector $\tilde{\boldsymbol{\epsilon}}(s)$ corrects for bias (underestimation of $\eta$ and overestimation of $\sigma$) as an extension of Banerjee et al. 2008, and $\mathcal{C}^T(\theta) \mathcal{C}^{*-1}(\theta)$ relates $w$ at desired point locations to the value of the latent GP at the knot locations, where $\mathcal{C}^T(\theta)$ is an $n \times m$ matrix that gets multiplied by the inverse of the $m \times m$ covariance matrix for the latent spatial process at the knots.
For a complete and general description see Finley et al. 2009, but here is the jist of the univariate model:

$$\boldsymbol{Y} \sim \text{Poisson}(\boldsymbol{X} \beta + \mathcal{C}^T(\theta) \mathcal{C}^{*-1}(\theta)\boldsymbol{w^*} + \tilde{\boldsymbol{\epsilon}})$$

$$\boldsymbol{w^*} \sim GP(0, \mathcal{C}^*(\theta))$$

$$\tilde{\boldsymbol{\epsilon}}(s) \sim \boldsymbol{N}(0, C(s, s) - \mathcal{c}^T(\theta) \mathcal{C}^{*-1}(\theta))\mathcal{c}(\theta)$$

with priors for $\eta$, $\sigma$, and $\phi$ completing a Bayesian specification.
In Stan syntax:


```
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
```

This approach seems to scale well for larger datasets relative to a full rank GP model. https://github.com/mbjoseph/gpp-speed-test
Comparing the number of effective samples per unit time for the two approaches across a range of sample sizes, it seems like that for larger datasets, the GPP executes more quickly and is more efficient ([code on GitHub](https://github.com/mbjoseph/gpp-speed-test)).



![](README_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Below are results from this small simulation study, where `n` is the number of data points in the sample, `elapsed_t` is the elapsed time in seconds, `n_eff` is the number of approximately independent samples for the log probability `lp__`, `Rhat` is a convergence diagnostic, and `efficiency` is `n_eff / elapsed_t`.


|   n| elapsed_t|    n_eff|     Rhat| efficiency|Model   |
|---:|---------:|--------:|--------:|----------:|:-------|
|  19|    21.062| 361.4927| 1.003165| 17.1632668|GPP     |
|  31|    31.883| 275.1871| 1.004865|  8.6311554|GPP     |
|  49|    55.072| 266.6816| 1.003334|  4.8424178|GPP     |
|  77|    79.444| 201.0303| 1.008528|  2.5304653|GPP     |
| 121|   126.924| 245.0213| 1.008804|  1.9304572|GPP     |
| 191|   242.799| 279.5732| 1.002820|  1.1514596|GPP     |
| 299|   449.470| 209.2244| 1.014566|  0.4654913|GPP     |
|  19|     4.347| 103.0432| 1.023866| 23.7044337|Full GP |
|  31|     8.869| 218.5975| 1.001799| 24.6473694|Full GP |
|  49|    21.174| 221.2929| 1.015297| 10.4511626|Full GP |
|  77|    66.045| 150.1677| 1.014973|  2.2737185|Full GP |
| 121|   296.052| 196.3052| 1.001702|  0.6630767|Full GP |
| 191|  1034.353| 228.9396| 1.008379|  0.2213361|Full GP |
| 299|  3459.850| 170.3018| 1.010420|  0.0492223|Full GP |

## More reading

[Finley, Andrew O., et al. "Improving the performance of predictive process modeling for large datasets." Computational statistics & data analysis 53.8 (2009): 2873-2884.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2743161/)

[Banerjee, Sudipto, et al. "Gaussian predictive process models for large spatial data sets." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.4 (2008): 825-848.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2741335/)

[Datta, Abhirup, et al. "Hierarchical nearest-neighbor Gaussian process models for large geostatistical datasets." Journal of the American Statistical Association. Accepted (2015).](http://arxiv.org/pdf/1406.7343.pdf)

[Guhaniyogi, Rajarshi, et al. "Adaptive Gaussian predictive process models for large spatial datasets." Environmetrics 22.8 (2011): 997-1007.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268014/)
