#7.5 Imputation: 
"The file interexp.dat contains data from an experiment that
was interrupted before all the data could be gathered. Of interest was the
difference in reaction times of experimental subjects when they were given
stimulus A versus stimulus B. Each subject is tested under one of the two
stimuli on their first day of participation in the study, and is tested under
the other stimulus at some later date. Unfortunately the experiment was
interrupted before it was finished, leaving the researchers with 26 subjects
with both A and B responses, 15 subjects with only A responses and 17
subjects with only B responses."

library(MASS)  # For multivariate normal functions

"a) Calculate empirical estimates of θA, θB, ρ, σ^2A, σ^2B from the data using
the commands mean , cor and var . Use all the A responses to get
ˆθA and ˆσ^2A, and use all the B responses to get ˆθB and ˆσ^2B. Use only
the complete data cases to get ˆρ."

# Load the data from GitHub
data_url <- "https://raw.githubusercontent.com/dstrick17/dstrickMA578/refs/heads/main/CASMA578/HW4/interexp.txt"
data <- read.table(data_url, header = TRUE, na.strings = "NA")

# Calculate empirical estimates
# a) Using all A responses for θA and σ^2_A
# na.rm means to only use data with non-missing values in that column
theta_A <- mean(data$yA, na.rm = TRUE)
sigma2_A <- var(data$yA, na.rm = TRUE)

# Using all B responses for θB and σ^2_B
theta_B <- mean(data$yB, na.rm = TRUE)
sigma2_B <- var(data$yB, na.rm = TRUE)

# Using complete cases (both yA and yB are non-NA) for ρ\
# Make data frame of only rows with data in both a and 
complete_cases <- data[complete.cases(data), ]
rho <- cor(complete_cases$yA, complete_cases$yB)
# Display results
cat("θA (Mean of A responses):", theta_A, "\n")
cat("θB (Mean of B responses):", theta_B, "\n")
cat("σ^2_A (Variance of A responses):", sigma2_A, "\n")
cat("σ^2_B (Variance of B responses):", sigma2_B, "\n")
cat("ρ (Correlation between A and B responses for complete cases):", rho, "\n")


#####################################################################
"b) For each person i with only an A response, impute a B response as
yˆi,B = (ˆθB + (yi,A −ˆθA)ˆρ * (σ^2B/σ^2A)^(1/2))
For each person i with only a B response, impute an A response as
yˆi,A = ˆθA + (yi,B − ˆθB)ˆρ* (σ^2A/σ^2B)^(1/2))
You now have two “observations” for each individual. Do a paired
sample t-test and obtain a 95% confidence interval for θA − θB."

# Impute missing B responses for those with only A responses
A_only <- is.na(data$yB) & !is.na(data$yA)
data$yB[A_only] <- theta_B + (data$yA[A_only] - theta_A) * rho * sqrt(sigma2_B / sigma2_A)

# Impute missing A responses for those with only B responses
B_only <- is.na(data$yA) & !is.na(data$yB)
data$yA[B_only] <- theta_A + (data$yB[B_only] - theta_B) * rho * sqrt(sigma2_A / sigma2_B)

# Paired sample t-test on the imputed data
t_test_result <- t.test(data$yA, data$yB, paired = TRUE)

# Print the t-test results and confidence interval for θA - θB
cat("Paired t-test results:\n")
print(t_test_result)
cat("\n95% Confidence interval for θA - θB:", t_test_result$conf.int, "\n")


##################################################################################
"c) Using either Jeffreys’ prior or a unit information prior distribution for
the parameters, implement a Gibbs sampler that approximates the
joint distribution of the parameters and the missing data. Compute
a posterior mean for θA − θB as well as a 95% posterior confidence
interval for θA − θB. Compare these results with the results from b)
and discuss."

# Set up prior parameters for the multivariate normal approach
mu_prior <- c(theta_A, theta_B)  # Prior mean vector for θ
L0 <- matrix(c(sigma2_A, rho * sqrt(sigma2_A * sigma2_B),
               rho * sqrt(sigma2_A * sigma2_B), sigma2_B), ncol = 2)  # Prior covariance matrix

nu0 <- 4  # Degrees of freedom for Wishart prior (p+2 where p=2 for bivariate)
S0 <- L0  # Prior scale matrix for the Wishart distribution

# Initialize variables
num_iter <- 10000
theta_samples <- matrix(0, nrow = num_iter, ncol = 2)  # Store θA and θB as columns in a matrix
sigma_samples <- array(0, dim = c(2, 2, num_iter))  # Store covariance matrices in 3D array

# Starting values
theta <- c(theta_A, theta_B)
Sigma <- L0  # Initial covariance matrix

# Ensure complete_cases is a numeric matrix
complete_cases <- as.matrix(complete_cases)

# Gibbs sampler loop
for (i in 1:num_iter) {
  
  # Step 1: Update θ (multivariate normal)
  n_complete <- nrow(complete_cases)  # Calculate the number of complete cases
  ybar <- colMeans(complete_cases)    # Mean of the complete cases
  
  # Posterior mean and covariance for θ
  L_n <- solve(solve(L0) + n_complete * solve(Sigma))
  mu_n <- L_n %*% (solve(L0) %*% mu_prior + n_complete * solve(Sigma) %*% ybar)
  
  # Sample from multivariate normal
  theta <- mvrnorm(1, mu = mu_n, Sigma = L_n)
  
  # Step 2: Update Sigma (inverse-Wishart using rWishart)
  # Center complete_cases by subtracting theta from each row
  centered_data <- sweep(complete_cases, 2, theta, "-")
  S_n <- S0 + t(centered_data) %*% centered_data
  
  # Sample from the inverse-Wishart by inverting the output of rWishart
  Sigma <- solve(rWishart(1, nu0 + n_complete, S_n)[,,1])  # Extract the first matrix
  
  # Step 3: Update missing data based on current θ and Sigma estimates
  A_only <- is.na(data$yB) & !is.na(data$yA)
  B_only <- is.na(data$yA) & !is.na(data$yB)
  
  # Impute missing B responses for those with only A responses
  data$yB[A_only] <- theta[2] + (data$yA[A_only] - theta[1]) * (Sigma[2, 1] / Sigma[1, 1])
  
  # Impute missing A responses for those with only B responses
  data$yA[B_only] <- theta[1] + (data$yB[B_only] - theta[2]) * (Sigma[1, 2] / Sigma[2, 2])
  
  # Store samples
  theta_samples[i, ] <- theta
  sigma_samples[,,i] <- Sigma
}



# Posterior summary
posterior_diff <- theta_samples[,1] - theta_samples[,2]
posterior_mean_diff <- mean(posterior_diff)
posterior_ci_diff <- quantile(posterior_diff, probs = c(0.025, 0.975))

# Print results
cat("Posterior mean for θA - θB:", posterior_mean_diff, "\n")
cat("95% Posterior confidence interval for θA - θB:", posterior_ci_diff, "\n")