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
parameters {
  // Core autoregressive model
  vector[N_train_years * N_sites] y;  // latent richness values
  real<lower=0,upper=1> beta_autoreg; // If β≥1, variance is not finite

  // regression model for the "anchor" (mean for mean-reversion)
  real alpha;
  vector[N_sites] alpha_site;
  vector[N_env] beta_env;

  // observation model
  vector[N_observers] alpha_observer;
  real<lower=0> sigma_observer;
  real<lower=0> sigma_error;

  // variance partitioning
  real<lower=0> sigma_total_latent;
  real<lower=0,upper=1> proportion_within;
}
transformed parameters {
  real sigma_within;
  real sigma_site;
  real sigma_autoreg;

  // Partition the total variance into within and between site variation
  sigma_within = sqrt(sigma_total_latent^2 * proportion_within);
  sigma_site = sqrt(sigma_total_latent^2 * (1 - proportion_within));

  // Determine how much to jump from year to year based on amount of
  // autocorrelation and the long-term variance of the AR(1) process.
  // See https://onlinecourses.science.psu.edu/stat510/node/60
  sigma_autoreg = sqrt(sigma_within^2 * (1 - beta_autoreg^2));
}
model {
  vector[num_elements(y)] env_effect;

  // priors on standard deviations
  sigma_error ~ gamma(2, 0.1);
  sigma_observer ~ gamma(2, 0.1);
  sigma_total_latent ~ gamma(2, 0.1);

  // Proportion of latent variation for intra-site variation
  proportion_within ~ beta(2, 2);

  // priors for regression model
  alpha ~ normal(0, 2);
  beta_env ~ normal(0, 2);

  // Autocorrelation
  beta_autoreg ~ beta(2, 2);

  // random effects
  alpha_site ~ normal(0, sigma_site);
  alpha_observer ~ normal(0, sigma_observer);

  if (include_environment) {
    env_effect = env * beta_env;
  } else {
    env_effect = rep_vector(0.0, num_elements(env_effect));
  }

  for (i in 1:N_sites) {
    int start;
    int end;
    vector[N_train_years] alpha_autoreg;
    vector[N_train_years] anchor;

    end = i * N_train_years;
    start = end - N_train_years + 1;

    // linear mixed model for long-term site-level means, plus perturbations
    // from the environment.
    anchor = alpha + alpha_site[i] + env_effect[start:end];

    // Calculate autoregressive alpha intercept from the long-term mean and
    // the autoregressive slope.
    // See https://onlinecourses.science.psu.edu/stat510/node/60
    alpha_autoreg = anchor * (1 - beta_autoreg);

    // First data point is drawn using long-term mean and sd
    y[start] ~ normal(anchor[1], sigma_within);

    // Subsequent data points depend on the year before: a + b*y
    y[(start+1):end] ~ normal(
            alpha_autoreg[2:N_train_years] + beta_autoreg * y[start:(end-1)],
            sigma_autoreg);
  }

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
      // Grab the final observation from this site
      previous_y = y[future_site_index[i] * N_train_years];
    } else{
      // Grab the value of y that was predicted for last year
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

    // Add the observer's bias and the observation noise
    future_observed[i] = normal_rng(future_y[i] + observer_effect, sigma_error);
  }
}
