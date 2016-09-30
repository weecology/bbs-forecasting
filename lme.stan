data {
  int N;
  int N_sites;
  int N_observers;
  vector[N] y;
  int site_index[N];
  int observer_index[N];
}
parameters {
  real intercept;
  real<lower=0> sigma;
  real<lower=0> sigma_site;
  real<lower=0> sigma_observer;
  vector[N_sites] raw_alpha_site;
  vector[N_observers] raw_alpha_observer;
}
transformed parameters {
  vector[N] mu;
  vector[N_sites] alpha_site;
  vector[N_observers] alpha_observer;
  
  // "Matt trick"
  alpha_site = raw_alpha_site * sigma_site;
  alpha_observer = raw_alpha_observer * sigma_observer;
  
  mu = intercept + alpha_site[site_index] + alpha_observer[observer_index];
}
model {
  intercept ~ normal(0, 5);
  
  sigma ~ exponential(1);
  sigma_site ~ exponential(1);
  sigma_observer ~ exponential(1);
  
  // "Matt trick"
  raw_alpha_site ~ normal(0, 1);
  raw_alpha_observer ~ normal(0, 1);
  
  y ~ normal(mu, sigma);
}
