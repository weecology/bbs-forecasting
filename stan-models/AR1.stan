data {
  // Sizes of vectors
  int<lower=0> N_observations; // Total non-NA observations in training set
  int<lower=0> N_train_years;
  int<lower=0> N_test_years;
  int<lower=0> N_sites;
  int<lower=0> N_observers;
  int<lower=0> N_env; // columns of the `env` matrix

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
  real<lower=0> sigma_diffusion;      // square root of diffusion rate
  real<lower=0> sigma_within;         // asymptotic within-site sd, I think

  // linear regression model for the "anchor" (mean for mean-reversion)
  real alpha;
  vector[N_sites] alpha_site_raw;
  vector[N_env] beta_env;
  real<lower=0> sigma_site;

  // observation model
  vector[N_observers] alpha_observer;
  real<lower=0> sigma_observer;
  real<lower=0> sigma_error;
}
transformed parameters {
  cov_matrix[N_train_years] Sigma;
  matrix[N_train_years, N_train_years] L_Sigma;

  for (i in 1:N_train_years) {
    for (j in 1:N_train_years) {
      // Covariance for a Weiner process ("standard Brownian motion"")
      Sigma[i,j] = min(i, j) * sigma_diffusion;
      if (i == j) {
        // Add a negligible amount to the diagonal to make the matrix more
        // invertible
        Sigma[i,j] = Sigma[i,j] + 1E-8;
      }
    }
  }

  // We're going to use this covariance matrix once per site, so it's helpful
  // to do the hard work associated with factorizing it here instead of inside
  // the loop.
  L_Sigma = cholesky_decompose(Sigma);
}
model {
  vector[num_elements(y)] env_effect;

  // priors on standard deviations
  sigma_site ~ gamma(2, 0.1);
  sigma_diffusion ~ gamma(2, 0.1);
  sigma_within ~ gamma(2, 0.1);
  sigma_error ~ gamma(2, 0.1);
  sigma_observer ~ gamma(2, 0.1);

  // priors for regression model on the environment
  alpha ~ normal(0, 2);
  beta_env ~ normal(0, 2);

  // random effects
  alpha_observer ~ normal(0, sigma_observer);
  // site-level effects weren't mixing well, so I'm trying the non-centering
  // approach described in the Stan manual (the "Matt trick"). The raw value
  // gets multiplied by sigma_site, which is equivalent to using sigma_site
  // as the sd here instead of 1.
  alpha_site_raw ~ normal(0, 1);

  if (include_environment) {
    env_effect = env * beta_env;
  } else {
    env_effect = rep_vector(0.0, num_elements(env_effect));
  }

  for (i in 1:N_sites) {
    int start; // where to find this site's first year
    int end;   // where to find this site's last year
    vector[N_train_years] y_start; // site's latent value in year 1
    vector[N_train_years] anchor; // Mean for mean reversion
    vector[N_train_years] alpha_autoreg; // intercept in linear autoregression
    real alpha_site;

    alpha_site = sigma_site * alpha_site_raw[i];

    end = i * N_train_years;
    start = end - N_train_years + 1;

    // linear mixed model for long-term site-level means, plus perturbations
    // from the environment.
    anchor = alpha + alpha_site + env_effect[start:end];

    // standard brownian motion (infinite asymptotic variance).
    // The mean for this component is just wherever the site started (y[start]).
    y_start = rep_vector(y[start], N_train_years);
    y[start:end] ~ multi_normal_cholesky(y_start, L_Sigma);
    // Mean-reversion to the anchor point
    y[start:end] ~ normal(anchor, sigma_within);
  }

  // Observation error; observed value is centered on y + observer effects
  scaled_richness ~ normal(
    y[observation_index] + alpha_observer[observer_index],
    sigma_error
  );
}
// generated quantities {
//   vector[N_sites * N_test_years] future_anchor;
//   vector[N_sites * N_test_years] future_alpha_autoreg;
//   vector[N_sites * N_test_years] future_y;
//   vector[N_sites * N_test_years] future_observed;
//
//   future_anchor = alpha + sigma_site * alpha_site_raw[future_site_index];
//   if (include_environment){
//     future_anchor = future_anchor + future_env * beta_env;
//   }
//   future_alpha_autoreg = future_anchor * (1 - beta_autoreg);
//
//   for (i in 1:num_elements(future_anchor)) {
//     real future_mu;
//     real previous_y;
//     real observer_effect;
//
//     if (i % N_test_years == 1) {
//       // Grab the final y value from this site's training years
//       previous_y = y[future_site_index[i] * N_train_years];
//     } else{
//       // Grab the value of y that was predicted for last year
//       previous_y = future_y[i - 1];
//     }
//
//     future_mu = future_alpha_autoreg[i] + beta_autoreg * previous_y;
//     future_y[i] = normal_rng(future_mu, sigma_autoreg);
//
//     // Add the effects of known observers; if the observer wasn't in the
//     // training set, then draw a random value for their observer effect.
//     if(future_observer_index[i] != 0){
//       observer_effect = alpha_observer[future_observer_index[i]];
//     } else{
//       observer_effect = normal_rng(0, sigma_observer);
//     }
//
//     // Add the observer's bias and the observation noise
//     future_observed[i] = normal_rng(future_y[i] + observer_effect, sigma_error);
//   }
// }
