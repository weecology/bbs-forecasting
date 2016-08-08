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

    // Autoregressive //
    y[(start+1):stop] ~ student_t(nu, beta * y[start:(stop-1)], sigma);
  }

  // Observation noise //
  observed ~ normal(y[index], nugget);
}
// generated quantities {
//   vector[N2] future_y[N_sites];
//   vector[N2] future_observed[N_sites];
//
//   for (i in 1:N_sites) {
//     for (j in 1:N2) {
//       if (j == 1) {
//         // First prediction depends on final year of y
//         future_y[i][1] = student_t_rng(nu, beta * y[i][N1], sigma);
//       } else {
//         // Subsequent predictions depend on preceding year //
//         future_y[i][j] = student_t_rng(nu, beta * future_y[i][j-1], sigma);
//       }
//
//       // Observation noise //
//       future_observed[i][j] = normal_rng(future_y[i][j], nugget);
//     }
//   }
// }
