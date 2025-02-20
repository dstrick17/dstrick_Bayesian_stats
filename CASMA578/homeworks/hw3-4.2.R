# Set the parameters for the posterior distributions of theta_A and theta_B
########### Exercise 3.3
# Tumor counts: A cancer laboratory is estimating the rate of tumorigenesis
# in two strains of mice, A and B. They have tumor count data for 10 mice
# in strain A and 13 mice in strain B. Type A mice have been well studied,
# and information from other laboratories suggests that type A mice have
# tumor counts that are approximately Poisson-distributed with a mean of
# 12. Tumor count rates for type B mice are unknown, but type B mice are
# related to type A mice. The observed tumor counts for the two populations
# are
# yA = (12, 9, 12, 14, 13, 13, 15, 8, 15, 6);
# yB = (11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7).
# a) Find the posterior distributions, means, variances and 95% quantilebased confidence intervals for θA and θB, assuming a Poisson sampling
# distribution for each group and the following prior distribution:
#   θA ∼ gamma(120,10), θB ∼ gamma(12,1), p(θA, θB) = p(θA)×p(θB). 4.2 Tumor count comparisons: Reconsider the tumor count data in Exercise
# 3.3:



############### 4.2 Tumor count comparisons: Reconsider the tumor count data in Exercise
# 3.3:
#   a) For the prior distribution given in part a) of that exercise, obtain
# Pr(θB < θA|yA, yB) via Monte Carlo sampling.

# b) For a range of values of n0, obtain Pr(θB < θA|yA, yB) for θA ∼
# gamma(120, 10) and θB ∼ gamma(12×n0, n0). Describe how sensitive
# the conclusions about the event {θB < θA} are to the prior distribution
# on θB.

# c) Repeat parts a) and b), replacing the event {θB < θA} with the event
# {Y˜B < Y˜A}, where Y˜A and Y˜B are samples from the posterior predictive distribution.
# ############################################################### Part A
# Set the parameters for the posterior distributions of theta_A and theta_B
shape_A <- 237  # Gamma shape parameter for theta_A
rate_A <- 20    # Gamma rate parameter for theta_A
shape_B <- 125  # Gamma shape parameter for theta_B
rate_B <- 14    # Gamma rate parameter for theta_B

# Number of Monte Carlo samples
n_samples <- 100000

# Generate samples from the posterior distributions
theta_A_samples <- rgamma(n_samples, shape = shape_A, rate = rate_A)
theta_B_samples <- rgamma(n_samples, shape = shape_B, rate = rate_B)

# Estimate the probability that theta_B < theta_A
# This function creates a vector of TRUE/FALSE values where element is true if 
# theta B is less than theta a for that sample in montecarlo simulation
# Taking the mean of the boolean vector gives the proportion of true values
prob_theta_B_less_than_theta_A <- mean(theta_B_samples < theta_A_samples)

# Calculate means and standard deviations of the posterior samples
mean_theta_A <- mean(theta_A_samples)
mean_theta_B <- mean(theta_B_samples)
sd_theta_A <- sd(theta_A_samples)
sd_theta_B <- sd(theta_B_samples)

# Output the results
cat("P(θB < θA | yA, yB) ≈", prob_theta_B_less_than_theta_A, "\n")
cat("Posterior mean of θA:", mean_theta_A, "\n")
cat("Posterior standard deviation of θA:", sd_theta_A, "\n")
cat("Posterior mean of θB:", mean_theta_B, "\n")
cat("Posterior standard deviation of θB:", sd_theta_B, "\n")

############################################################# Part B

# Range of n0 values to explore
n0_values <- seq(0.5, 10, by = 0.5)

# Initialize a vector to store probabilities for each n0
prob_theta_B_less_than_theta_A_list <- numeric(length(n0_values))

# Loop over each value of n0
for (i in seq_along(n0_values)) {
  
  # Set the shape and rate for theta_B based on the current value of n0
  n0 <- n0_values[i]
  shape_B <- 12 * n0  # Gamma shape parameter for theta_B
  rate_B <- n0        # Gamma rate parameter for theta_B
  
  # Generate samples from the posterior distributions
  theta_A_samples <- rgamma(n_samples, shape = shape_A, rate = rate_A)
  theta_B_samples <- rgamma(n_samples, shape = shape_B, rate = rate_B)
  
  # Calculate the probability that theta_B < theta_A
  prob_theta_B_less_than_theta_A <- mean(theta_B_samples < theta_A_samples)
  
  # Store the result
  prob_theta_B_less_than_theta_A_list[i] <- prob_theta_B_less_than_theta_A
}

# Output the results
cat("n0 values:\n", n0_values, "\n")
cat("Probabilities of θB < θA:\n", prob_theta_B_less_than_theta_A_list, "\n")

# Plot the results
plot(n0_values, prob_theta_B_less_than_theta_A_list, type = "b", pch = 19, col = "blue",
     xlab = "n0", ylab = "Pr(θB < θA | yA, yB)",
     main = "Sensitivity of Pr(θB < θA) to the Prior on θB",
     ylim = c(0, 1))

######################################################################## Part C

# Generate samples from the posterior predictive distributions
# Poisson samples using theta_A and theta_B from their respective posteriors
Y_A_samples <- rnbinom(n_samples, size = shape_A, mu = theta_A_samples)
Y_B_samples <- rnbinom(n_samples, size = shape_B, mu = theta_B_samples)

# Estimate the probability that Y_B < Y_A
prob_Y_B_less_than_Y_A <- mean(Y_B_samples < Y_A_samples)

# Calculate means and standard deviations of the posterior predictive samples
mean_Y_A <- mean(Y_A_samples)
mean_Y_B <- mean(Y_B_samples)
sd_Y_A <- sd(Y_A_samples)
sd_Y_B <- sd(Y_B_samples)

# Output the results
cat("P(Y_B < Y_A | yA, yB) ≈", prob_Y_B_less_than_Y_A, "\n")
cat("Posterior predictive mean of Y_A:", mean_Y_A, "\n")
cat("Posterior predictive standard deviation of Y_A:", sd_Y_A, "\n")
cat("Posterior predictive mean of Y_B:", mean_Y_B, "\n")
cat("Posterior predictive standard deviation of Y_B:", sd_Y_B, "\n")

##################################################################### Part C for Range of N0

# Range of n0 values to explore
n0_values <- seq(0.5, 10, by = 0.5)

# Initialize a vector to store probabilities for each n0
prob_Y_B_less_than_Y_A_list <- numeric(length(n0_values))

# Loop over each value of n0
for (i in seq_along(n0_values)) {
  
  # Set the shape and rate for theta_B based on the current value of n0
  n0 <- n0_values[i]
  shape_B <- 12 * n0  # Gamma shape parameter for theta_B
  rate_B <- n0        # Gamma rate parameter for theta_B
  
  # Generate samples from the posterior distributions
  theta_A_samples <- rgamma(n_samples, shape = shape_A, rate = rate_A)
  theta_B_samples <- rgamma(n_samples, shape = shape_B, rate = rate_B)

  
  # Calculate the probability that Y_B < Y_A
  prob_Y_B_less_than_Y_A <- mean(Y_B_samples < Y_A_samples)
  
  # Store the result
  prob_Y_B_less_than_Y_A_list[i] <- prob_Y_B_less_than_Y_A
}

# Output the results
cat("n0 values:\n", n0_values, "\n")
cat("Probabilities of Y_B < Y_A:\n", prob_Y_B_less_than_Y_A_list, "\n")

# Plot the results
plot(n0_values, prob_Y_B_less_than_Y_A_list, type = "b", pch = 19, col = "blue",
     xlab = "n0", ylab = "Pr(Y_B < Y_A | yA, yB)",
     main = "Sensitivity of Pr(Y_B < Y_A) to the Prior on θB",
     ylim = c(0, 1))
