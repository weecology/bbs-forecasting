data {
  // Sizes of vectors
  int<lower=0> N_observations; // Total non-NA observations
  int<lower=0> N_train_years;  // Number of training years
  int<lower=0> N_sites;
  int<lower=0> N_observers;

  // Pointers to information about each observation
  int<lower=0> observation_index[N_observations];
  int<lower=0> observer_index[N_observations];
  int<lower=0> site_index[N_sites * N_train_years];

  // Pointers to specific kinds of years within sites
  int<lower=0> which_first[N_sites];
  int<lower=0> which_non_first[N_sites * (N_train_years - 1)];

  // Predictor variables
  vector[N_sites * N_train_years] ndvi_sum;

  // response variable
  real scaled_richness[N_observations];
}
transformed data {
  int<lower=0> which_non_last[N_sites * (N_train_years - 1)];

  // Apparently I need a loop for this
  for (i in 1:size(which_non_first)) {
    which_non_last[i] = which_non_first[i] - 1;
  }
}
parameters {
  // Core autoregressive model
  vector[N_train_years * N_sites] y;
  real<lower=0> sigma_autoreg;
  real beta_autoreg;

  // regression model for the "anchor" (mean for mean-reversion)
  real alpha;
  real<lower=0> sigma_site;
  vector[N_sites] alpha_site;
  real beta_ndvi;

  // observation model
  real<lower=0> sigma_observer;
  real<lower=0> sigma_error;
  vector[N_observers] alpha_observer;
}
transformed parameters {
  vector[N_train_years * N_sites] anchor;
  vector[N_sites * N_train_years] alpha_autoreg;
  vector[(N_train_years - 1) * N_sites] mu_non_first;

  // Anchor (mean for mean-reversion) for each site/year combination.
  // This is basically a linear mixed model, with one level per site.
  anchor = alpha + alpha_site[site_index] + beta_ndvi * ndvi_sum;

  // convert from "anchor & beta" to more standard "alpha & beta"
  alpha_autoreg = anchor * (1 - beta_autoreg);

  // autoregression: determine expected value of y for next year based on
  // current value of y
  mu_non_first = alpha_autoreg[which_non_last] +
      beta_autoreg * y[which_non_last];
}
model {
  // priors on standard deviations
  sigma_autoreg  ~ gamma(2, 0.01);
  sigma_error ~ gamma(2, 0.01);
  sigma_site ~ gamma(2, 0.01);
  sigma_observer ~ gamma(2, 0.01);

  // Prior on regression intercept
  alpha ~ normal(0, 10);

  // Prior on autoregressive beta (probably between 0 and 1)
  beta_autoreg ~ normal(0.5, 0.5);

  // random effects
  alpha_site ~ normal(0, sigma_site);
  alpha_observer ~ normal(0, sigma_observer);

  // Weak prior on y1 for when no observations were taken in year 1
  //    (sd==2 is a very weak prior when the data is scaled to have sd==1)
  y[which_first] ~ normal(anchor[which_first], 2);

  // Autoregression model for non-first years: inferred values are centered
  // on mu
  y[which_non_first] ~ normal(mu_non_first, sigma_autoreg);

  // Observation error; observed value is centered on y + observer effects
  scaled_richness ~ normal(
    y[observation_index] + alpha_observer[observer_index],
    sigma_error
  );
}

