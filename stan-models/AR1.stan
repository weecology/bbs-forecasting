data {
  // Sizes of vectors
  int<lower=0> N_observations; // Total non-NA observations in training set
  int<lower=0> N_train_years;
  int<lower=0> N_test_years;
  int<lower=0> N_sites;
  int<lower=0> N_observers;
  int<lower=0> N_env;

  // Pointers to information about each observation
  int<lower=0> observation_index[N_observations];
  int<lower=0> observer_index[N_observations];
  int<lower=0> site_index[N_sites * N_train_years];

  // Pointers to specific kinds of years within sites
  int<lower=0> which_first[N_sites];
  int<lower=0> which_non_first[N_sites * (N_train_years - 1)];

  // Predictor variables
  matrix[N_sites * N_train_years, N_env] env;

  // response variable
  real scaled_richness[N_observations];

  // Flag to include environmental predictors
  int<lower=0,upper=1> include_environment;

  // From the future
  matrix[N_sites * N_test_years, N_env] future_env;
  int future_site_index[N_sites * N_test_years];
  int future_observer_index[N_sites * N_test_years];
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
  vector[N_env] beta_env;

  // observation model
  real<lower=0> sigma_observer;
  real<lower=0> sigma_error;
  vector[N_observers] alpha_observer;
}
transformed parameters {
  vector[N_train_years * N_sites] anchor;
  // Anchor (mean for mean-reversion) for each site/year combination.
  // This is basically a linear mixed model, with one level per site.
  anchor = alpha + alpha_site[site_index];
  if (include_environment){
    anchor = anchor + env * beta_env;
  }
}
model {
  vector[N_sites * N_train_years] alpha_autoreg;
  vector[(N_train_years - 1) * N_sites] mu_non_first;

  // priors on standard deviations
  sigma_autoreg  ~ gamma(2, 0.1);
  sigma_error ~ gamma(2, 0.1);
  sigma_site ~ gamma(2, 0.1);
  sigma_observer ~ gamma(2, 0.1);

  // priors for regression model
  alpha ~ normal(0, 10);
  beta_env ~ normal(0, 2);

  beta_autoreg ~ normal(0.5, 0.5); // (probably between 0 and 1)

  // random effects
  alpha_site ~ normal(0, sigma_site);
  alpha_observer ~ normal(0, sigma_observer);

  // convert from "anchor & beta" to more standard "alpha & beta"
  alpha_autoreg = anchor * (1 - beta_autoreg);

  // autoregression: determine expected value of y for next year based on
  // current value of y
  mu_non_first = alpha_autoreg[which_non_last] +
      beta_autoreg * y[which_non_last];

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
generated quantities {
  vector[N_sites * N_test_years] future_anchor;
  vector[N_sites * N_test_years] future_alpha_autoreg;
  vector[N_sites * N_test_years] future_y;
  vector[N_sites * N_test_years] future_observed;

  future_anchor = alpha + alpha_site[future_site_index];
  if (include_environment){
    future_anchor = future_anchor + future_env * beta_env;
  }
  future_alpha_autoreg = future_anchor * (1 - beta_autoreg);

  for (i in 1:num_elements(future_anchor)) {
    real future_mu;
    real previous_y;
    real observer_effect;

    if (i % N_test_years == 1) {
      previous_y = y[future_site_index[i] * N_train_years];
    } else{
      previous_y = future_y[i - 1];
    }

    future_mu = future_alpha_autoreg[i] + beta_autoreg * previous_y;
    future_y[i] = normal_rng(future_mu, sigma_autoreg);

    // Add the effects of known observers; if the observer wasn't in the
    // training set, then draw a random value for their observer effect.
    if(future_observer_index[i] != 0){
      observer_effect = alpha_observer[future_observer_index[i]];
    } else{
      observer_effect = normal_rng(0, sigma_observer);
    }
    future_observed[i] = normal_rng(future_y[i] + observer_effect, sigma_error);
  }
}
