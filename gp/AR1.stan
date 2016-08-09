data {
  int<lower=0> N_obs; // Total non-NA observations
  int<lower=0> N1;    // Number of training years
  int<lower=0> N2;    // Number of test years
  int<lower=0> N_sites;

  int<lower=0> index[N_obs];
  real observed[N_obs];
}
parameters {
  // scalars //
  real<lower=0> sigma;
  real<lower=0> nugget;
  real<lower=0> nu_minus_two; // see transformed parameters
  real beta;

  // Latent variables //
  vector[N1 * N_sites] y;
}
transformed parameters {
  // If nu is less than 2, variance is infinite
  real<lower=2> nu;
  nu = nu_minus_two + 2;
}
model {
  // priors //
  sigma  ~ gamma(3, 0.1);
  nugget ~ gamma(3, 0.1);
  nu_minus_two ~ gamma(3, 0.1);

  beta ~ normal(0, 1);

  // likelihood //
  for (i in 1:N_sites){
    int start;
    int stop;
    start = (i - 1) * N1 + 1;
    stop = i * N1;

    y[(start+1):stop] ~ student_t(nu, beta * y[start:(stop-1)], sigma);
  }

  observed ~ normal(y[index], nugget);
}
generated quantities {
  vector[N2 * N_sites] future_y;
  vector[N2 * N_sites] future_observed;

  for (i in 1:N_sites){
    int start;
    start = (i - 1) * N2 + 1;

    for(j in 1:N2){
      real mu_temp;
      int to_update;
      to_update = start + j - 1;

      if (j == 1){
        // First prediction is based on last value of future_y
        mu_temp = beta * y[i * N1];
      } else{
        // subsequent predictions are based on previous value of future_y
        mu_temp = beta * future_y[to_update - 1];
      }

      future_y[to_update] = student_t_rng(nu, mu_temp, sigma);
      future_observed[to_update] = normal_rng(future_y[to_update], nugget);
    }
  }
}
