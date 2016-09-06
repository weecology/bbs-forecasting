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
  // scalars //
  real<lower=0> sigma_autoreg;
  real<lower=0> sigma_observer;
  real<lower=0> sigma_site;
  real<lower=0> sigma_site_beta;
  real<lower=0> nugget;

  // Autoregressive model //
  vector[N_sites] site_alpha;
  vector[N_sites] site_beta;
  real beta;

  // Latent variables //
  vector[N_train_years * N_sites] y;
  vector[N_observers] observer_alpha;
}
transformed parameters {
  real mu_non_first[(N_train_years - 1) * N_sites];

  // For non-first data points, what was richness the preceding year? //
  for (i in 1:num_elements(mu_non_first)) {
    int n;
    // alpha + beta * y, where alpha and beta are for the current _site_, and
    // y is for the previous _year_
    n = site_index[which_non_first[i]];
    mu_non_first[i] = site_alpha[n] +
      (beta + site_beta[n]) * y[which_non_first[i] - 1];
  }
}
model {
  // priors on standard deviations //
  sigma_autoreg  ~ gamma(2, 0.01);
  nugget ~ gamma(2, 0.01);
  sigma_site ~ gamma(2, 0.01);
  sigma_site_beta ~ gamma(2, 0.01);
  sigma_observer ~ gamma(2, 0.01);

  // prior on autoregressive rate //
  beta ~ normal(0.5, 0.5); // probably between 0 and 1

  // random effects //
  site_beta ~ normal(0, sigma_site_beta);
  site_alpha ~ normal(0, sigma_site);
  observer_alpha ~ normal(0, sigma_observer);

  // Weak prior on y1 for when no observations were taken in year 1 //
  y[which_first] ~ normal(site_alpha[site_index[which_first]], 2);

  // Autoregression, with mean reverting to zero. //
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
