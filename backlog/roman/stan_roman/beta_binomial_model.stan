// estimate proportion
data {
  int n;
  int y;
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  theta ~ beta(8, 160);
  y ~ binomial(n, theta);
}

generated quantities {
  real<lower=0, upper=1> prior_theta;
  prior_theta = beta_rng(8, 160);
}








