data {
  // Sizes of vectors //
  int<lower=0> N_observations; // Total non-NA observations
  int<lower=0> N_train_years;  // Number of training years
  int<lower=0> N_sites;
  int<lower=0> N_observers;

  // Pointers to information about each observation //
  int<lower=0> observation_index[N_observations];
  int<lower=0> observer_index[N_observations];
  int<lower=0> site_index[N_sites * N_train_years];

  // Pointers to specific years within sites//
  int<lower=0> which_first[N_sites];
  int<lower=0> which_non_first[N_sites * (N_train_years - 1)];

  // response variable (Z-scaled) //
  real scaled_richness[N_observations];
}
parameters {
  // Core autoregressive process model //
  real<lower=0> sigma_autoreg;
  real<lower=0> nugget;
  vector[N_train_years * N_sites] y;

  // site-level model //
  //    First vector is site-level mean
  //    Second vector is memory (autoregressive "beta").
  corr_matrix[2] site_cor;
  vector[2] site_coefs[N_sites];
  vector[2] mu_site;
  vector<lower=0>[2] sigma_site;

  // observation model //
  real<lower=0> sigma_observer;
  vector[N_observers] observer_alpha;
}
transformed parameters {
  vector[N_sites] site_means;
  vector[N_sites] site_alphas;
  vector[N_sites] site_betas;

  real mu_non_first[(N_train_years - 1) * N_sites];

  for (i in 1:N_sites) {
    site_means[i] = site_coefs[i, 1];
    site_betas[i] = site_coefs[i, 2];
    site_alphas[i] = site_means[i] * (1 - site_betas[i]);
  }

  // For non-first data points, what was richness the preceding year? //
  for (i in 1:num_elements(mu_non_first)) {
    int n;
    // alpha + beta * y, where alpha and beta are for the _current site_,
    // and y is for the _previous year_
    n = site_index[which_non_first[i]];
    mu_non_first[i] = site_alphas[n] +
      site_betas[n] * y[which_non_first[i] - 1];
  }
}
model {
  // priors on standard deviations //
  sigma_autoreg  ~ gamma(2, 0.01);
  nugget ~ gamma(2, 0.01);
  sigma_site ~ gamma(2, 0.01);
  sigma_observer ~ gamma(2, 0.01);

  // prior on global mean richness //
  mu_site[1] ~ normal(0, 0.1);     // must be close to observed global mean
  // prior on autoregressive rate
  mu_site[2] ~ normal(0.5, 0.5);   // probably between 0 and 1

  // random effects //
  site_coefs ~ multi_normal(mu_site, quad_form_diag(site_cor, sigma_site));
  observer_alpha ~ normal(0, sigma_observer);

  // Weak prior on y1 for when no observations were taken in year 1 //
  // (sd==2 is a weak prior when the data is scaled to have sd==1)  //
  y[which_first] ~ normal(site_means[site_index[which_first]], 2);

  // Autoregression //
  y[which_non_first] ~ normal(mu_non_first, sigma_autoreg);

  // Observation error //
  scaled_richness ~ normal(
    y[observation_index] + observer_alpha[observer_index],
    nugget
  );
}

