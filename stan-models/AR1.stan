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

  // response variable
  real scaled_richness[N_observations];
}
parameters {
  // Core autoregressive model
  vector[N_train_years * N_sites] y;
  real<lower=0> sigma_autoreg;

  // site-level model //
  //    First vector is site-level mean
  //    Second vector is memory (autoregressive "beta").
  corr_matrix[2] site_cor;
  vector[2] site_coefs[N_sites];
  vector[2] mu_site;
  vector<lower=0>[2] sigma_site;

  // observation model //
  real<lower=0> sigma_observer;
  real<lower=0> sigma_error;
  vector[N_observers] observer_alpha;
}
transformed parameters {
  vector[N_sites] site_means;
  vector[N_sites] site_alphas;
  vector[N_sites] site_betas;
  real mu_non_first[(N_train_years - 1) * N_sites];

  // Extract and simplify coefficients.
  for (i in 1:N_sites) {
    site_means[i] = site_coefs[i, 1];
    site_betas[i] = site_coefs[i, 2];
    // long-term expected value of an AR1 process
    // is alpha / (1-beta)
    site_alphas[i] = site_means[i] * (1 - site_betas[i]);
  }

  // Define the autoregressive means based on previous year's y.
  for (i in 1:num_elements(mu_non_first)) {
    int site_num;
    site_num = site_index[which_non_first[i]];
    // mu = alpha + beta * y, where alpha and beta are for the _current site_,
    // and y is for the _previous year_
    mu_non_first[i] = site_alphas[site_num] +
      site_betas[site_num] * y[which_non_first[i] - 1];
  }
}
model {
  // priors on standard deviations
  sigma_autoreg  ~ gamma(2, 0.01);
  sigma_error ~ gamma(2, 0.01);
  sigma_site ~ gamma(2, 0.01);
  sigma_observer ~ gamma(2, 0.01);

  // prior on global mean richness (probably close to 0 after scalilng)
  mu_site[1] ~ normal(0, 0.1);
  // prior on autoregressive rate (probably between 0 and 1)
  mu_site[2] ~ normal(0.5, 0.5);

  // random effects
  site_coefs ~ multi_normal(mu_site, quad_form_diag(site_cor, sigma_site));
  observer_alpha ~ normal(0, sigma_observer);

  // Weak prior on y1 for when no observations were taken in year 1
  //    (sd==2 is a very weak prior when the data is scaled to have sd==1)
  y[which_first] ~ normal(site_means[site_index[which_first]], 2);

  // Autoregression for non-first years; mu defined in "transformed parameters"
  y[which_non_first] ~ normal(mu_non_first, sigma_autoreg);

  // Observation error
  scaled_richness ~ normal(
    y[observation_index] + observer_alpha[observer_index],
    sigma_error
  );
}

