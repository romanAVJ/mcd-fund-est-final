###############################################################################
# raw script for final exam
# roman_avj 2023/11/30
###############################################################################
library(tidyverse)

# ex. 1 -------------------------------------------------------------------
#### 1.1 (us proportion) ####
# params
PROP_HIGH_SCHOOL <- 0.26 # claim
num_high_school <- 69 # observed number of high school graduates
num_sample <- 310 # sample size
prop_investigator <- num_high_school / num_sample # sample proportion

# does the data support the claim? #
# H0: p = 0.26
# H1: p < 0.26

# test statistic
se_prop <- sqrt(prop_investigator * (1 - prop_investigator) / num_sample)
z_stat <- (prop_investigator - PROP_HIGH_SCHOOL) / se_prop # wald test
# p-value
pval <- pnorm(z_stat)

# conclusion: the p_val is near 0.05 so it is not very clear if the data
# supports the claim. The p-value is not less than 0.05 so classically we do not reject the null hypothesis.
# Maybe using bayesian methods we can get a better answer.

#### 1.2 (mendell's peas) ####
# params
P_PEAS <- c(9/16, 3/16, 3/16, 1/16) # claim
observed_peas <- c(315, 101, 108, 32) # observed number of peas

# does the data support the claim? #
# H0: p = P_PEAS
# H1: p != P_PEAS
# Use likelihood ratio test
create_multin_loglike <- function(x){
  # x is a vector of observed counts
  # returns a function that calculates the log likelihood of a multinomial
  # distribution with parameters
  function(p){
    # p is a vector of parameters
    # returns the log likelihood of a multinomial distribution with parameters p
    # and observed counts x
    sum(x * log(p))
  }
}

# get log likelihoods #
# null hyp
loglike_null <- create_multin_loglike(observed_peas)(P_PEAS)
# alt hyp
p_maxlike <- observed_peas / sum(observed_peas)
loglike_alt <- create_multin_loglike(observed_peas)(p_maxlike)

# test statistic #
# likelihood ratio test
chi_stat = 2 * (loglike_alt - loglike_null)
# p-value
pval <- pchisq(chi_stat, df = 3, lower.tail = FALSE)

# conclusion: the pval is 0.92, no evidence to reject the null hypothesis.
# Very high p-value, no necessity to do a parametric bootstrap

#### 1.3 (wald test for poisson) ####
#### 1.3.1 (theoretic statistic)
wald_test <- function(x, lambda){
  # x is the observed number of events
  # lambda is the claim
  # returns a list with the test statistic and the p-value
  # H0: lambda = lambda
  # H1: lambda != lambda
  # test statistic
  lambda_estimate <- mean(x)
  se_lambda <- sqrt(lambda_estimate / length(x))
  z_stat <- (lambda_estimate - lambda) / se_lambda
  # p-value
  pval <- 2 * pnorm(-abs(z_stat))
  return(pval)
}

#### 1.3.2 (empirical statistic)
# params
LAMBDA0 <- 1
N_SAMPLE <- 20
ALPHA <- 0.05
N_BOOT <- 1000
SEED <- 8

# function to simulate data and use test
simulate_wald <- function(lambda0, n_sample, n_boot, seed){
  # lambda0 is the claim
  # n_sample is the sample size
  # n_boot is the number of bootstrap samples
  # seed is the seed for the random number generator
  # returns a list with the test statistic and the p-value
  # generate data
  set.seed(seed)
  # bootstrap
  pval_boot <- replicate(
    n=n_boot,
    expr=wald_test(
      x=rpois(n_sample, lambda0),
      lambda=lambda0
      )
  )
  # return
  return(pval_boot)
}

# generate data
pval_wald_poiss <- simulate_wald(
  lambda0=LAMBDA0,
  n_sample=N_SAMPLE,
  n_boot=N_BOOT,
  seed=SEED
)

# look number of times we reject the null hypothesis
mean(pval_wald_poiss < ALPHA)

# conclusion: the number of times we reject the null hypothesis is really near to 0.05,
# so the type I error is near 0.05


# bootstrap and bayesian inference ----------------------------------------
#### 2.1 (parametric bootstrap of a normal) ####
# read data
load("data/x.RData")
# rename x to x_obs
x_obs <- x
# delete x
rm(x)

#### 2.1.1 log like of a normal with known mean
create_normal_loglike <- function(x){
  # x is a vector of observed values
  # returns a function that calculates the log likelihood of a normal
  # distribution with known mean
  function(sigma2){
    # sigma is the standard deviation
    # returns the log likelihood of a normal distribution with known mean
    # and standard deviation sigma
    -length(x) * log(sigma2) - sum((x - 0)^2) / (sigma2)
  }
}

# get log likelihood
loglike_xobs <- create_normal_loglike(x_obs)

#### 2.1.2 get max likelihood
# get max like by optimization
sigma2_maxlike <- optim(
  par=1,
  fn=loglike_xobs,
  method="Brent",
  control=list(fnscale=-1),
  lower=0,
  upper=1000
)$par


# plot log likelihood
sigma2 <- seq(1, 200, 1)
plot(sigma2, loglike_xobs(sigma2), type="l")
# ad vline in max likelihood
abline(v=sigma2_maxlike, col="red")
# comment: my the graph, it seems that the maximum likelihood is not very
# clear where the maximum is. We are inferring for sigma^2, not sigma.

#### 2.1.3 parametric bootstrap
simulate_sigma2 <- function(sigma2, n_obs, n_boot, seed=8){
  # n_boot is the number of bootstrap samples
  # sigma2 is the claim
  # simulate data
  set.seed(seed)
  # return bootstrap
  sigma2_boot <- replicate(
    n=n_boot,
    expr=var(rnorm(n_obs, mean=0, sd=sqrt(sigma2)))
  )
}
# simulate bootstrap
sigma2_boot <- simulate_sigma2(
  sigma2=sigma2_maxlike,
  n_obs=length(x_obs),
  n_boot=10000,
  seed=8
)
# get standard error of the bootstrap
se_sigma2_boot <- sd(sigma2_boot)

# plot bootstrap and add vertical line in claim
hist(sigma2_boot, breaks=20)
abline(v=sigma2_maxlike, col="red")

#### 2.2 bayesian analysis ####
#### 2.2.1 (prior: inverse gamma)
# params
# i am not completly sure of the typical value, then alpha must be small
ALPHA <- 1.5
beta <- ALPHA * sigma2_maxlike

# plot prior for tau
taus <- seq(1, 1000, 1)
plot(taus, dgamma(taus, shape=ALPHA, scale=beta), type="l")

# plot prior for sigma2
dinv_gamma <- function(sigma2, alpha, beta){
  # sigma2 is the standard deviation
  # alpha is the shape
  # beta is the scale
  # returns the density of the inverse gamma distribution
  alpha * beta^alpha / gamma(alpha) * sigma2^(-alpha - 1) * exp(-beta / sigma2)
}
sigmas2 <- seq(1, 1000, 1)
plot(sigmas2, dinv_gamma(taus, alpha=ALPHA, beta=beta), type="l")


#### 2.2.2 (analyticaly get the posterior)
# get posterior (manually did it)

#### 2.2.3 (simulate from posterior)
# simulate from posterior
simulate_sigma2_posterior <- function(x, alpha, beta, seed=8, n_sim=1000){
  # x is the observed data
  # alpha is the shape
  # beta is the scale
  # seed is the seed for the random number generator

  # update alpha and beta
  alpha_posterior <- alpha + length(x) / 2
  beta_posterior <- beta + sum(x^2) / 2

  # simulate from posterior
  set.seed(seed)
  sigma2_posterior <- 1 / rgamma(n=n_sim, shape=alpha_posterior, rate=beta_posterior)
  return(sigma2_posterior)
}

# simulate
sigma2_posterior <- simulate_sigma2_posterior(
  x=x_obs,
  alpha=ALPHA,
  beta=beta,
  seed=8,
  n_sim=10000
)

# se
se_sigma2_posterior <- sd(sigma2_posterior)


# graph posterior and add vertical line in claim
hist(sigma2_posterior, breaks=20)
abline(v=sigma2_maxlike, col="red")

#### 2.2.4 (compare bootstrap and bayesian)
# compare distributions
df_compar <- tibble(
    sigma2_boot=sigma2_boot,
    sigma2_posterior=sigma2_posterior
  ) |>
  pivot_longer(
    cols=c(sigma2_boot, sigma2_posterior),
    names_to="method",
    values_to="sigma2"
  ) |>
  mutate(
    method=fct_recode(method, "boot"="sigma2_boot", "bayes"="sigma2_posterior")
  )

# plot densities
ggplot(df_compar, aes(x=sigma2, fill=method)) +
  geom_density(alpha=0.5) +
  theme_minimal()

# plot ecdf
ggplot(df_compar, aes(x=sigma2, color=method)) +
  stat_ecdf() +
  theme_minimal()

# qq plots between bootstrap and bayesian
x_boot <- df_compar |> filter(method == "boot") |> pull(sigma2) |> sort()
x_posterior <- df_compar |> filter(method == "bayes") |> pull(sigma2) |> sort()
qqplot(x_boot, x_posterior)
# add diagonal
abline(a=0, b=1, col="red")

#### 2.3 (inference over function of parameters) ####
#### 2.3.1
# inference over t = log(sigma)
# what is the maximum likelihood of t?
# t_mle = log(sigma2_mle) = log(sd(x_obs))
t_mle <- log(sqrt(sigma2_maxlike))

# use parametric bootstrap to get the standard error of t_mle
simulate_t <- function(sigma, n_obs, n_boot, seed=8){
  # n_obs is the number of observations
  # n_boot is the number of bootstrap samples
  # seed is the seed for the random number generator

  # simulate bootstrap
  set.seed(seed)
  t_mle_boot <- replicate(
    n=n_boot,
    expr=log(sd(rnorm(n_obs, mean=0, sd=sigma)))
  )
  return(t_mle_boot)
}
# simulate t
t_sim <- simulate_t(
  sigma=sqrt(sigma2_maxlike),
  n_obs=length(x_obs),
  n_boot=10000,
  seed=8
)
# get 95% confidence interval
alpha_ci <- 0.05
t_freq_ci <- quantile(t_sim, probs=c(alpha_ci/2, 1-alpha_ci/2))

# graph t
hist(t_sim, breaks=20)
abline(v=t_mle, col="red")
# add 95% confidence interval
abline(v=t_freq_ci, col="red", lty=2)

#### 2.3.2
# bayesian inference over t
# simulate from posterior
simulate_t_posterior <- function(x, alpha, beta, seed=8, n_sim=1000){
  # x is the observed data
  # alpha is the shape
  # beta is the scale
  # seed is the seed for the random number generator

  # update alpha and beta
  alpha_posterior <- alpha + length(x) / 2
  beta_posterior <- beta + sum(x^2) / 2

  # simulate from posterior
  set.seed(seed)
  sigma2_posterior <- 1 / rgamma(n=n_sim, shape=alpha_posterior, rate=beta_posterior)
  t_posterior <- log(sqrt(sigma2_posterior))
  return(t_posterior)
}

# simulate
t_posterior <- simulate_t_posterior(
  x=x_obs,
  alpha=ALPHA,
  beta=beta,
  seed=8,
  n_sim=10000
)

# get 95% confidence interval
alpha_ci <- 0.05
t_bayes_ci <- quantile(t_posterior, probs=c(alpha_ci/2, 1-alpha_ci/2))

# plot posterior and add vertical line in claim
hist(t_posterior, breaks=20)
abline(v=t_mle, col="red")
# add 95% confidence interval
abline(v=t_bayes_ci, col="red", lty=2)

#### compare bootstrap and bayesian
# compare distributions
df_compar <- tibble(
    t_boot=t_sim,
    t_posterior=t_posterior
  ) |>
  pivot_longer(
    cols=c(t_boot, t_posterior),
    names_to="method",
    values_to="t"
  ) |>
  mutate(
    method=fct_recode(method, "boot"="t_boot", "bayes"="t_posterior")
  )

# plot densities
ggplot(df_compar, aes(x=t, fill=method)) +
  geom_density(alpha=0.5) +
  theme_minimal()

# plot ecdf
ggplot(df_compar, aes(x=t, color=method)) +
  stat_ecdf() +
  theme_minimal()

# they look quite similar

# ex 3. bayesian regularization -------------------------------------------
#### 3.1: eda ####
# read an dta file
df_polls <- foreign::read.dta("data/pew_research_center_june_elect_wknd_data.dta")
glimpse(df_polls)
# read election file
df_elections <- read_csv("data/2008ElectionResult.csv")
glimpse(df_elections)

#### 3.1.1: proportion of polls very liberal
## get tables
# prop as very liberals
table_polls <- df_polls |>
  mutate(
    is_very_liberal = (ideo == "very liberal")
  ) |>
  group_by(state) |>
  summarise(
    total_vliberal = sum(is_very_liberal, na.rm=TRUE),
    prop_vliberal = mean(is_very_liberal, na.rm=TRUE),
    n_polls = sum(!is.na(is_very_liberal))
  ) |>
  arrange(state) |>
  filter(
    !(state %in% c("alaska", "hawaii", "washington dc")) # remove these states
  )

# prop of votes for obama
table_obama <- df_elections |>
  mutate(state = tolower(state)) |>
  inner_join(table_polls, by=c("state")) |>
  select(state, prop_vliberal, vote_Obama_pct)
table_obama

## plot
# graph scatter between n_polls and prop_vliberal, add the name of the state
table_polls |>
  mutate(state_abr = str_sub(state, 1, 2)) |>
  ggplot(aes(x=n_polls, y=prop_vliberal)) +
    geom_point() +
    geom_text(aes(label=state_abr), hjust=0, vjust=0) +
    theme_minimal()


# graph scatter between prop_vliberal and vote_Obama_pct, add the name of the state
table_obama |>
  mutate(state_abr = str_sub(state, 1, 2)) |>
  ggplot(aes(x=prop_vliberal, y=vote_Obama_pct)) +
    geom_point() +
    geom_text(aes(label=state), hjust=0, vjust=0) +
    theme_minimal()

#### 3.2: bayesian inference ####
#### 3.2.1 estimate proportion using bayes stats ####
# params
A <- 8
B <- 160

# simulate posterior of a proportion
simulate_prop_posterior <- function(n_success, n_obs, n_sim=1000, A=2, B=2, seed=8){
  # n_obs is the number of observations
  # n_sim is the number of simulations
  # A is the number of successes
  # B is the number of failures
  # seed is the seed for the random number generator

  # update A and B
  A_posterior <- A + n_success
  B_posterior <- B + (n_obs - n_success)

  # simulate from posterior
  set.seed(seed)
  prop_posterior <- rbeta(n=n_sim, shape1=A_posterior, shape2=B_posterior)
  return(prop_posterior)
}

# simulate posterior for each state
df_obama_posterior <- table_polls |>
  # simulate and save them in a list and then expand
  group_by(state) |>
  mutate(
    prop_posterior = list(simulate_prop_posterior(
      n_success=total_vliberal,
      n_obs=n_polls,
      n_sim=10000,
      A=A,
      B=B,
      seed=8
    ))
  ) |>
  ungroup() |>
  select(state, prop_posterior) |>
  unnest(prop_posterior)

# get mean, std and 95% confidence interval for each state
table_stats_obama <- df_obama_posterior |>
  group_by(state) |>
  summarise(
    mean = mean(prop_posterior),
    std = sd(prop_posterior),
    ci_lower = quantile(prop_posterior, probs=0.025),
    ci_upper = quantile(prop_posterior, probs=0.975)
  ) |>
  arrange(state)

#### graph
# ploy distribution for each state
df_obama_posterior |>
  ggplot(aes(x=prop_posterior)) +
  geom_histogram(bins=50) +
  geom_vline(
    data=table_stats_obama,
    aes(xintercept=mean),
    color="red"
  ) +
  geom_vline(
    data=table_stats_obama,
    aes(xintercept=ci_lower),
    color="red",
    linetype="dashed"
  ) +
  geom_vline(
    data=table_stats_obama,
    aes(xintercept=ci_upper),
    color="red",
    linetype="dashed"
  ) +
  facet_wrap(~state, ncol=5) +
  theme_minimal()

#### 3.2.2 use STAN for idaho & virginia ####
# get stan model
model_stan <- file.path("stan_roman", "beta_binomial_model.stan") |> cmdstanr::cmdstan_model()
model_stan

## get data for idaho
data_idaho <- table_polls |>
  filter(state == 'idaho')
data_idaho

# fit model
fit_stan <- model_stan$sample(
  data = list(
    n = data_idaho$n_polls,
    y = data_idaho$total_vliberal
  ),
  chains = 8,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 10000,
  refresh = 500,
  seed = 8
)

# diagnose
fit_stan$cmdstan_diagnose()
fit_stan$summary()
fit_stan$draws(c("theta", "prior_theta")) |> 
  posterior::as_draws_df() |> 
  ggplot(aes(.iteration, theta)) +
  geom_line() +
  facet_wrap(~.chain, ncol = 1)

fit_stan$draws(c("theta", "prior_theta")) |> 
  posterior::as_draws_df() |> 
  mutate(.chain = factor(.chain)) |>
  ggplot(aes(.iteration, theta, color = .chain)) +
  geom_line(alpha = 0.5)

## get data from virginia
data_virginia <- table_polls |>
  filter(state == 'virginia')
data_virginia

# fit model
fit_stan <- model_stan$sample(
  data = list(
    n = data_virginia$n_polls,
    y = data_virginia$total_vliberal
  ),
  chains = 8,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 10000,
  refresh = 500,
  seed = 8
)

# diagnose
fit_stan$cmdstan_diagnose()
fit_stan$summary()
fit_stan$draws(c("theta", "prior_theta")) |>
  posterior::as_draws_df() |> 
  ggplot(aes(.iteration, theta)) +
  geom_line() +
  facet_wrap(~.chain, ncol = 1)

fit_stan$draws(c("theta", "prior_theta")) |>
  posterior::as_draws_df() |> 
  mutate(.chain = factor(.chain)) |>
  ggplot(aes(.iteration, theta, color = .chain)) +
  geom_line(alpha = 0.5)


#### 3.2.3 repeate 3.2.1 with point estimate ####
# join data
table_obama_posterior  <- table_polls |>
  inner_join(table_stats_obama, by = c("state")) |>
  mutate(
    prop_vliberal = mean
  ) 

## plot
# graph scatter between n_polls and prop_vliberal, using the point estimate from table_stats_obama
table_obama_posterior |> 
  mutate(method = "bayesian") |> 
  select(state, total_vliberal, prop_vliberal, n_polls, method) |>
  add_row(
    table_polls |> mutate(method = "frequentist") |> select(state, total_vliberal, prop_vliberal, n_polls, method)
  ) |> 
  mutate(state_abr = str_sub(state, 1, 2)) |>
  ggplot(aes(x=n_polls, y=prop_vliberal, color=method)) +
  geom_point() +
  geom_text(aes(label=state_abr), hjust=0, vjust=0) +
  theme_minimal()
# note: the proportion looks more controlled with the bayesian method

# graph scatter between prop_vliberal and vote_Obama_pct, using the point estimate from table_stats_obama
df_elections |>
  mutate(state = tolower(state)) |>
  inner_join(table_obama_posterior, by=c("state")) |>
  select(state, prop_vliberal, vote_Obama_pct) |>
  mutate(state_abr = str_sub(state, 1, 2)) |>
  ggplot(aes(x=prop_vliberal, y=vote_Obama_pct)) +
    geom_point() +
    geom_text(aes(label=state), hjust=0, vjust=0) +
    theme_minimal()
