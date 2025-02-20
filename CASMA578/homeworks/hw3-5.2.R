# 5.2 Sensitivity analysis: Thirty-two students in a science classroom were
# randomly assigned to one of two study methods, A and B, so that
# nA = nB = 16 students were assigned to each method. After several
# weeks of study, students were examined on the course material with an
# exam designed to give an average score of 75 with a standard deviation of
# 10. The scores for the two groups are summarized by {y¯A = 75.2, sA = 7.3}
# and {y¯B = 77.5, sb = 8.1}. Consider independent, conjugate normal prior
# distributions for each of θA and θB, with µ0 = 75 and σ2 = 100 for
# both groups. For each (κ0, ν0) ∈ {(1,1),(2,2),(4,4),(8,8),(16,16),(32,32)}
# (or more values), obtain Pr(θA < θB|yA, yB) via Monte Carlo sampling.
# Plot this probability as a function of (κ0 = ν0). Describe how you might
# use this plot to convey the evidence that θA < θB to people of a variety
# of prior opinions.

############### Define  Data and parameters for group A and B
n_A <- 16
n_B <- 16

ybar_A <- 75.2 # Sample mean for A
s_A <- 7.3 # Sample standard deviation for A
ybar_B <- 77.5
s_B <- 8.1

# Prior parameters
mu_0 <- 75
sigma0_sq <- 100
kappa_values <- c(1, 2, 4, 8, 16, 32) # sensitivity values for kappa0 and nu0


################ Find posterior mean and variance
posterior_mean <- function(kappa0, ybar, n, mu0) {
  (kappa0 * mu0 + n * ybar) / (kappa0 + n)
}

posterior_variance <- function(kappa0, sigma0_sq, n) {
  sigma0_sq / (kappa0 + n)
}


############### Use MonteCarlo sampling
set.seed(123) # for reproducibility
monte_carlo_pr_theta_A_lt_theta_B <- function(kappa0, n_A, ybar_A, s_A, n_B, ybar_B, s_B, mu_0, sigma0_sq, n_samples = 10000) {
  # Compute posterior mean and variance for group A and B
  mean_A <- posterior_mean(kappa0, ybar_A, n_A, mu_0)
  var_A <- posterior_variance(kappa0, sigma0_sq, n_A)
  
  mean_B <- posterior_mean(kappa0, ybar_B, n_B, mu_0)
  var_B <- posterior_variance(kappa0, sigma0_sq, n_B)
  
  # Generate Monte Carlo samples
  # rnorm is a built in function in R that generates random samples from a normal dist.
  samples_A <- rnorm(n_samples, mean_A, sqrt(var_A))
  samples_B <- rnorm(n_samples, mean_B, sqrt(var_B))
  
  # Estimate Pr(theta_A < theta_B) - creates True False vector recording how many A is less than B
  mean(samples_A < samples_B)
}

# Run the Monte Carlo sampling for each kappa0 value
results <- sapply(kappa_values, function(kappa0) {
  monte_carlo_pr_theta_A_lt_theta_B(kappa0, n_A, ybar_A, s_A, n_B, ybar_B, s_B, mu_0, sigma0_sq)
})

# Display results
data.frame(kappa_values, results)



############### Plot Results
library(ggplot2)
df <- data.frame(kappa_values, results)
ggplot(df, aes(x = kappa_values, y = results)) +
  geom_line() +
  geom_point() +
  labs(title = "Sensitivity Analysis: Pr(θA < θB)",
       x = "κ0 = ν0",
       y = "Pr(θA < θB)") +
  theme_minimal()

