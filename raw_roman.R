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
plot(sigma, loglike_xobs(sigma2), type="l")
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
# get posterior






























