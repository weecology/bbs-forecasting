data {
  int N;
  int N_site;
  int N_observer;
  int site_index[N];
  int observer_index[N];
  real richness[N];
}
parameters {
  vector[N_site] site_effect;
  vector[N_observer] observer_effect;
  real intercept;
  real<lower=0> site_sigma;
  real<lower=0> observer_sigma;
  real<lower=0> sigma;
}
model {
  // priors
  intercept ~ normal(mean(richness), 5 * sd(richness));
  
  site_sigma ~ normal(0, sd(richness));
  observer_sigma ~ normal(0, sd(richness));
  sigma ~ normal(0, sd(richness));
  
  // Latent variables
  site_effect ~ normal(0, site_sigma);
  observer_effect ~ normal(0, observer_sigma);
  
  // observation model
  richness ~ normal(
    intercept + site_effect[site_index] + observer_effect[observer_index], 
    sigma
  );
}
generated quantities {
  // Model prediction for how many species the average observer would find
  vector[N] expected_richness;
  expected_richness = to_vector(richness) - observer_effect[observer_index];
}
