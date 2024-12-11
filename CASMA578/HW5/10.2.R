# 10.2 Nesting success: Younger male sparrows may or may not nest during a
# mating season, perhaps depending on their physical characteristics. Researchers 
# have recorded the nesting success of 43 young male sparrows
# of the same age, as well as their wingspan, and the data appear in the
# file msparrownest.dat. Let Yi be the binary indicator that sparrow i
# successfully nests, and let xi denote their wingspan. Our model for Yi
# is logit Pr(Yi = 1|α, β, xi) = α + βxi, where the logit function is given by
# logit θ = log[θ/(1 − θ)].

# a) Write out the joint sampling distribution Qn
# i=1 p(yi|α, β, xi) and simplify as much as possible.

# b) Formulate a prior probability distribution over α and β by considering the 
# range of Pr(Y = 1|α, β, x) as x ranges over 10 to 15, the
# approximate range of the observed wingspans.


# 
library(tidyverse)
library(coda)

# Data preparation
sparrow_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrick_Bayesian_stats/refs/heads/main/CASMA578/HW5/msparrownest.txt", header = TRUE)
colnames(sparrow_data) <- c("y", "x")  # Rename columns for clarity

# Y is binary variable indicating if the sparrow has successfully nested
# X is the wing length

# Inspect data
summary(sparrow_data)

# Define log-likelihood function
log_likelihood <- function(beta, gamma, data) {
  y <- data$y
  x <- data$x
  eta <- beta + gamma * x
  p <- 1 / (1 + exp(-eta))
  sum(y * log(p) + (1 - y) * log(1 - p))
}

# Define log-prior function
log_prior <- function(beta, gamma) {
  dnorm(beta, mean = 0, sd = 10, log = TRUE) + dnorm(gamma, mean = 0, sd = 10, log = TRUE)
}

# Define log-posterior function
log_posterior <- function(beta, gamma, data) {
  log_likelihood(beta, gamma, data) + log_prior(beta, gamma)
}

##########################################################################
#########################################################################
# c) Implement a Metropolis algorithm that approximates p(α, β|y, x).
# Adjust the proposal distribution to achieve a reasonable acceptance
# rate, and run the algorithm long enough so that the effective sample
# size is at least 1,000 for each parameter.
# Metropolis sampler

metropolis_sampler <- function(data, n_iter = 10000, init = c(0, 0), proposal_sd = c(0.1, 0.1)) {
  # Initialize parameters
  beta_current <- init[1]
  gamma_current <- init[2]
  samples <- matrix(NA, nrow = n_iter, ncol = 2)
  colnames(samples) <- c("beta", "gamma")
  
  # Run the sampler
  for (i in 1:n_iter) {
    # Propose new values
    beta_proposed <- rnorm(1, mean = beta_current, sd = proposal_sd[1])
    gamma_proposed <- rnorm(1, mean = gamma_current, sd = proposal_sd[2])
    
    # Compute acceptance ratio
    log_r <- log_posterior(beta_proposed, gamma_proposed, data) - 
      log_posterior(beta_current, gamma_current, data)
    r <- exp(log_r)
    
    # Accept or reject
    if (runif(1) < r) {
      beta_current <- beta_proposed
      gamma_current <- gamma_proposed
    }
    
    # Store current samples
    samples[i, ] <- c(beta_current, gamma_current)
  }
  
  as.data.frame(samples)
}

# Run the sampler
set.seed(42)
samples <- metropolis_sampler(sparrow_data, n_iter = 10000, proposal_sd = c(0.1, 0.1))


# d) Compare the posterior densities of α and β to their prior densities.

# e) Using output from the Metropolis algorithm, come up with a way to
# make a confidence band for the following function fαβ(x) of wingspan:
# fαβ(x) = e α+βx1 + eα+βx ,
# where α and β are the parameters in your sampling model. Make a
# plot of such a band.

# Convert samples to mcmc object
library(coda)
mcmc_samples <- as.mcmc(samples)

# Check convergence diagnostics
summary(mcmc_samples)
plot(mcmc_samples)

# Effective sample size
effectiveSize(mcmc_samples)

# Summarize posterior samples
summary(samples)

# Visualize posterior densities
library(ggplot2)

ggplot(samples, aes(x = beta)) +
  geom_density(fill = "blue", alpha = 0.3) +
  labs(title = "Posterior Density of Beta", x = "Beta", y = "Density")

ggplot(samples, aes(x = gamma)) +
  geom_density(fill = "green", alpha = 0.3) +
  labs(title = "Posterior Density of Gamma", x = "Gamma", y = "Density")

# Define logistic function
f <- function(beta, gamma, x) {
  exp(beta + gamma * x) / (1 + exp(beta + gamma * x))
}

# Generate confidence bands
x_range <- seq(10, 15, length.out = 100)
ci_band <- apply(samples, 1, function(params) f(params[1], params[2], x_range))
ci_band <- apply(ci_band, 1, quantile, probs = c(0.025, 0.975))

# Plot confidence bands
plot(x_range, apply(ci_band, 2, mean), type = "l", col = "blue", ylim = range(ci_band),
     xlab = "Wingspan (x)", ylab = "Probability of Nesting", main = "Confidence Band for Logistic Function")
lines(x_range, ci_band[1, ], col = "red", lty = 2)
lines(x_range, ci_band[2, ], col = "red", lty = 2)
legend("topright", legend = c("Mean", "95% CI"), col = c("blue", "red"), lty = c(1, 2))

