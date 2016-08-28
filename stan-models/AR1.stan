data {
  int<lower=0> N_observations; // Total non-NA observations
  int<lower=0> N1;    // Number of training years
  int<lower=0> N2;    // Number of test years
  int<lower=0> N_sites;
  int<lower=0> N_observers;

  int<lower=0> index[N_observations];
  int<lower=0> observer_index[N_observations];
  real observed[N_observations];
}
parameters {
  // scalars //
  real<lower=0> sigma;
  real<lower=0> nugget;
  real<lower=0> nu;
  vector[N_sites] alpha;
  real beta;
  real<lower=0> observer_sigma;
  real<lower=0>sigma_alpha;

  // Latent variables //
  vector[N1 * N_sites] y;
  vector[N_observers] observer_alpha;
}
transformed parameters {
}
model {
  // priors //
  sigma  ~ gamma(3, 0.1);
  nugget ~ gamma(3, 0.1);
  observer_sigma ~ gamma(3, 0.1);
  sigma_alpha ~ gamma(3, 0.1);

  nu ~ gamma(3, 0.1);

  alpha ~ normal(0, sigma_alpha);
  observer_alpha ~ normal(0, observer_sigma);

  beta ~ normal(0.5, 0.5); // probably between 0 and 1

  // likelihood //
  for (i in 1:N_sites){
    int start;
    int stop;
    start = (i - 1) * N1 + 1;
    stop = i * N1;

    // Weak prior on y1 for when no observations were taken in year 1
    y[start] ~ student_t(nu, alpha[i], 2 * sd(observed));

    // Autoregressive component
    y[(start+1):stop] ~ student_t(
      nu,
      alpha[i] + beta * y[start:(stop-1)],
      sigma
    );
  }

  // Observation error
  observed ~ normal(y[index] + observer_alpha[observer_index], nugget);
}
generated quantities {
  vector[N2 * N_sites] future_y;
  //vector[N2 * N_sites] future_observed;

  for (i in 1:N_sites){
    int start;
    start = (i - 1) * N2 + 1;

    for(j in 1:N2){
      real mu_temp;
      int to_update;
      to_update = start + j - 1;

      if (j == 1){
        // First prediction is based on last value of future_y
        mu_temp = alpha[i] + beta * y[i * N1];
      } else{
        // subsequent predictions are based on previous value of future_y
        mu_temp = alpha[i] + beta * future_y[to_update - 1];
      }

      future_y[to_update] = student_t_rng(nu, mu_temp, sigma);
      //future_observed[to_update] = normal_rng(future_y[to_update], nugget);
    }
  }
}
