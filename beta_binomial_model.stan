// beta_binomial_model.stan
data {
  int<lower=0> n_success; // número de éxitos
  int<lower=0> n_trials;  // número de ensayos
  real<lower=0> alpha_prior; // parámetro alpha de la Beta a priori
  real<lower=0> beta_prior;  // parámetro beta de la Beta a priori
}
parameters {
  real<lower=0, upper=1> p; // proporción de éxitos
}
model {
  p ~ beta(alpha_prior, beta_prior); // prior Beta
  n_success ~ binomial(n_trials, p);  // likelihood
}
