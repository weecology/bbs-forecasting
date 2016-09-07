data {
  // Sizes of vectors //
  int<lower=0> N_observations; // Total non-NA observations
  int<lower=0> N_train_years;  // Number of training years
  int<lower=0> N_test_years;   // Number of test years
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
  real observed_richness[N_observations];
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
    // alpha + beta * y, where alpha and beta are for the _current site_, and
    // y is for the _previous year_
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
  // prior on autoregressive rate  //
  mu_site[2] ~ normal(0.5, 0.5);   // probably between 0 and 1

  // random effects //
  site_coefs ~ multi_normal(mu_site, quad_form_diag(site_cor, sigma_site));
  observer_alpha ~ normal(0, sigma_observer);

  // Weak prior on y1 for when no observations were taken in year 1 //
  y[which_first] ~ normal(site_means[site_index[which_first]], 2);

  // Autoregression //
  y[which_non_first] ~ normal(mu_non_first, sigma_autoreg);

  // Observation error
  observed_richness ~ normal(y[observation_index] + observer_alpha[observer_index], nugget);
}
// generated quantities {
//   vector[N_test_years * N_sites] future_y;
//   //vector[N_test_years * N_sites] future_observed;
//
//   for (i in 1:N_sites){
//     int start;
//     start = (i - 1) * N_test_years + 1;
//
//     for(j in 1:N_test_years){q
//       real mu_temp;
//       int to_update;
//       to_update = start + j - 1;
//
//       if (j == 1){
//         // First prediction is based on last value of future_y
//         mu_temp = site_trend[i] + beta * y[i * N_train_years];
//       } else{
//         // subsequent predictions are based on previous value of future_y
//         mu_temp = site_trend[i] + beta * future_y[to_update - 1];
//       }
//
//       // future_y[to_update] = student_t_rng(nu, mu_temp +, sigma_autoreg);
//       //future_observed[to_update] = normal_rng(future_y[to_update], nugget);
//     }
//   }
// }
