# 3.3 part A
# Prior is beta(2,8), find theta for teen recidivism

n <- 43 # number of teens released from prison
y <- 15 # Number of re-offenders going back to prison within 36 months

alpha_prior_A <- 2
beta_prior_A <- 8

alpha_posterior_A <- alpha_prior_A + y
beta_posterior_A <- beta_prior_A + n

posterior_mean_A <- alpha_posterior_A / (alpha_posterior_A + beta_posterior_A)

posterior_mode_A <- (beta_posterior_A - 1) / (alpha_posterior_A + beta_posterior_A - 2)

post_stdev_A <- sqrt((alpha_posterior_A * beta_posterior_A) / ((alpha_posterior_A + beta_posterior_A)^2 *(alpha_posterior_A + beta_posterior_A + 1)))

ci_A <- qbeta(c(0.025, 0.975), alpha_posterior_A, beta_posterior_A)

# Results
cat("Part (a) Results:\n")
cat("Posterior Mean:", posterior_mean_A, "\n")
cat("Posterior Mode:", posterior_mode_A, "\n")
cat("Posterior SD:", post_stdev_A, "\n")
cat("95% CI:", ci_A, "\n\n")

#############################
# 3.4 part B

# Prior is beta(8,2), find theta for teen recidivism

alpha_prior_B <- 8
beta_prior_B <- 2

alpha_posterior_B <- alpha_prior_B + y
beta_posterior_B <- beta_prior_B + n

posterior_mean_B <- alpha_posterior_B / (alpha_posterior_B + beta_posterior_B)

posterior_mode_B <- (beta_posterior_B - 1) / (alpha_posterior_B + beta_posterior_B - 2)

post_stdev_B <- sqrt((alpha_posterior_B * beta_posterior_B) / ((alpha_posterior_B + beta_posterior_B)^2 *(alpha_posterior_B + beta_posterior_B + 1)))

ci_B <- qbeta(c(0.025, 0.975), alpha_posterior_B, beta_posterior_B)

# Results
cat("Part (b) Results:\n")
cat("Posterior Mean:", posterior_mean_B, "\n")
cat("Posterior Mode:", posterior_mode_B, "\n")
cat("Posterior SD:", post_stdev_B, "\n")
cat("95% CI:", ci_B, "\n\n")

#############################################
# 3.4 part c

# Function to calculate density of mixtures
mixture_density <- function(x) {
  0.75 * dbeta(x, 2, 8) + 0.25 * dbeta(x, 8, 2)
}

# Create sequence of values between 0 and 1
x <- seq(0, 1, length.out = 1000)

# Calculate the density for the mixture and the individual beta distributions
mixture_densities <- mixture_density(x)
beta_A_densities <- dbeta(x, 2, 8)
beta_B_densities <- dbeta(x, 8, 2)

# Calculate max density for y-axis limits
max_density <- max(mixture_densities, beta_A_densities, beta_B_densities)

# Plot mixture distribution with adjusted y-axis limits
plot(x, mixture_densities, type = "l", col = "purple", lwd = 3,
     xlab = "Probability of re-offending", ylab = "Density",
     main = "Comparison of Prior Distributions",
     ylim = c(0, max_density * 1.1))  # Increase limit by 10%

# Add the beta(2,8) and beta(8,2) distributions to the plot
lines(x, beta_A_densities, col = "blue", lwd = 2)
lines(x, beta_B_densities, col = "red", lwd = 2)

# Add legend
legend("top", c("Mixture", "Beta(2,8)", "Beta(8,2)"),
       col = c("purple", "blue", "red"), lwd = 2)

########################################################
# 3.4 part d

# Likelihood function
likelihood <- function(theta) {
  dbinom(y, n, theta)
}

# Unnormalized posterior (prior * likelihood)
unnormalized_posterior <- function(theta) {
  mixture_density(theta) * likelihood(theta)
}

# Generate theta values
theta_values <- seq(0, 1, length.out = 1000)

# Calculate unnormalized posterior for each theta
posterior_values <- sapply(theta_values, unnormalized_posterior)

# Plot
plot(theta_values, posterior_values, type = "l", 
     xlab = "theta", ylab = "p(theta) × p(y|theta)", 
     main = "Unnormalized Posterior Distribution")

# Find the mode (maximum of unnormalized posterior)
mode_index <- which.max(posterior_values)
posterior_mode <- theta_values[mode_index]

cat("Approximate posterior mode:", posterior_mode, "\n")

# Add mode to the plot
abline(v = posterior_mode, col = "red", lty = 2)
text(posterior_mode, max(posterior_values), "Mode", pos = 4, col = "red")

# Modes for Beta(a,b) 
mode_a <- (2 + y - 1) / (2 + 8 + n - 2)  # Beta(2,8) prior
mode_b <- (8 + y - 1) / (8 + 2 + n - 2)  # Beta(8,2) prior

cat("Mode for part a (Beta(2,8) prior):", mode_a, "\n")
cat("Mode for part b (Beta(8,2) prior):", mode_b, "\n")

# Add these modes to the plot
abline(v = mode_a, col = "blue", lty = 2)
abline(v = mode_b, col = "green", lty = 2)
legend("topright", c("Mixture Posterior", "Beta(2,8) Posterior", "Beta(8,2) Posterior"),
       col = c("red", "blue", "green"), lty = 2)

