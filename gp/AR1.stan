data {
  int<lower=0> N1;
  int<lower=0> N_sites;
  vector[N1] y1[N_sites];
}
parameters {
  // scalars //
  real<lower=0> sigma;
  real<lower=0> nugget;
  real<lower=0> nu_minus_two; // see transformed parameters
  real beta;

  // Latent variables //
  vector[N1] y[N_sites];
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
    // Autoregressive //
    tail(y[i], N1 - 1) ~ student_t(nu, beta * head(y[i], N1 - 1), sigma);

    // Observation noise //
    y1[i] ~ normal(y[i], nugget);
  }
}
