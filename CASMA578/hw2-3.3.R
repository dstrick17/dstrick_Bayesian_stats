####################
# 3.3 part A
# p(θA, θB) = p(θA)×p(θB)
# Number of tumors in mice in sample A
yA = c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6)
nA <- length(yA)
sum_Ya <- sum(yA)

# Theta_A is gamma(120,10)
# Prior parameters for for θA
aA <- 120; bA <- 10

# Posterior parameters for θA
posterior_aA <- aA + sum_Ya # Shape parameter
posterior_bA <- bA + nA # Rate parameter
mean_θA <- posterior_aA / posterior_bA # Mean
variance_θA <- posterior_aA / (posterior_bA^2) # Variance
ci_θA <- qgamma(c(0.025, 0.975), shape = posterior_aA, rate = posterior_bA)# 95% Confidence Interval 


# Number of tumors in mice of sample B
yB <- c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7)
nB <- length(yB)
sum_yB <- sum(yB)

# Theta_B is gamma(12,1)
# Prior parameters for θB
aB <- 12; bB <- 1

# Posterior parameters for θB
posterior_aB <- aB + sum_yB # Shape
posterior_bB <- bB + nB # Rate
mean_θB <- posterior_aB / posterior_bB # Mean = shape/rate
variance_θB <- posterior_aB / (posterior_bB^2) # Variance
ci_θB <- qgamma(c(0.025, 0.975), shape = posterior_aB, rate = posterior_bB)# 95% Confidence Interval 

# Results
cat("Posterior distribution for θA:\n")
cat("Mean:", mean_θA, "\n","Variance:", variance_θA, "\n","95% CI", ci_θA, "\n")

cat("Posterior distribution for θB:\n")
cat("Mean:", mean_θB, "\n", "Variance:",variance_θB, "\n","95% CI",ci_θB, "\n")

##############################################
#3.3 part B


# Initialize parameters, a larger n0 places more weight on the prior
n0_values <- 1:50  # Values for n0

posterior_expectations <- sapply(n0_values, function(n0) {
  a_prior <- 12 * n0
  b_prior <- n0
  (a_prior + sum_yB) / (b_prior + nB)
})

plot(n0_values, posterior_expectations, type = "l", 
     xlab = "n0", ylab = "Posterior Expectation of θB",
     main = "Posterior Expectation of θB vs. Prior distribution (n0)")

# In order for the posterior expectation of θB to be close to that of θA,
# Prior parameters should reflect the belief that θB is close to θA and 
# The evidence (marginal likelihood) should not overwhelmingly contradict the prior belief

##############################################
# 3.3 part c

# If p(θA, θB) = p(θA)×p(θB), then population A and population B are 
# conditionally independent so knowledge about population A
# would not affect our beliefs about population B

#####################################################################
######################################################################



