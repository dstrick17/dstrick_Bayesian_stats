# 6.2 Mixture model:
"The file glucose.dat contains the plasma glucose concentration of 532 females 
 from a study on diabetes (see Exercise 7.6)."

# Load the glucose data from GitHub
glucose_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrickMA578/main/CASMA578/HW4/glucose.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)

# Check the first few rows to confirm the structure of the data
head(glucose_data)


"a) Make a histogram or kernel density estimate of the data. Describe
 how this empirical distribution deviates from the shape of a normal
 distribution"

# Rename the column 
colnames(glucose_data) <- "glucose"

# Create a histogram of the glucose concentrations
hist(glucose_data$glucose, 
     main = "Histogram of Glucose Concentrations", 
     xlab = "Glucose Concentration", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black")

"b) Consider the following mixture model for these data: For each study
participant there is an unobserved group membership variable Xi
which is equal to 1 or 2 with probability p and 1 − p. If Xi = 1
then Yi ∼ normal(θ1, σ1**U2), and if Xi = 2 then Yi ∼ normal(θ2, σ2**2). Let
p ∼ beta(a, b), θj ∼ normal(µ0,τ0)**2 and 1/σj ∼ gamma(ν0/2, ν0σ0**2/2
for both j = 1 and j = 2. Obtain the full conditional distributions of
(X1, . . . , Xn), p, θ1, θ2, σ1**2 and σ2**2"

# See answers on paper


"c) Setting a = b = 1, µ0 = 120, τ0^2= 200, σ0^2 = 1000 and ν0 = 10,
implement the Gibbs sampler for at least 10,000 iterations. Let
θ1^(s) = min{θ1^(s), θ2^(s)} and θ2^(s) = max{θ1^(s), θ2^(s)}. 
Compute and plot the autocorrelation functions of θ1^(s) and θ2^(s), 
as well as their effective sample sizes."

# Load necessary library
library(MCMCpack)  # For the inverse gamma function
library(coda)  # For autocorrelation and effective sample size calculations

# Load the data
glucose_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrickMA578/refs/heads/main/CASMA578/HW4/glucose.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)
y <- glucose_data$V1
n <- length(y)

# Prior settings
a <- 1; b <- 1
mu0 <- 120; tau0_sq <- 200
sigma0_sq <- 1000; nu0 <- 10

# Initial values for Gibbs sampler
theta1 <- mean(y) + rnorm(1)
theta2 <- mean(y) - rnorm(1)
sigma1_sq <- var(y)
sigma2_sq <- var(y)
p <- 0.5
X <- sample(1:2, n, replace = TRUE)

# Number of iterations
num_iter <- 10000

# Store samples
theta1_samples <- numeric(num_iter)
theta2_samples <- numeric(num_iter)
p_samples <- numeric(num_iter)
sigma1_sq_samples <- numeric(num_iter)
sigma2_sq_samples <- numeric(num_iter)
theta_min_samples <- numeric(num_iter)
theta_max_samples <- numeric(num_iter)

# Gibbs Sampler
for (s in 1:num_iter) {
  
  # Step 1: Sample each Xi
  prob1 <- p * dnorm(y, mean = theta1, sd = sqrt(sigma1_sq))
  prob2 <- (1 - p) * dnorm(y, mean = theta2, sd = sqrt(sigma2_sq))
  probs <- prob1 / (prob1 + prob2)
  X <- rbinom(n, size = 1, prob = probs) + 1
  
  # Step 2: Update p
  n1 <- sum(X == 1)
  n2 <- n - n1
  p <- rbeta(1, a + n1, b + n2)
  
  # Step 3: Update theta1 and theta2
  y1 <- y[X == 1]
  y2 <- y[X == 2]
  sigma1_post <- 1 / (1 / tau0_sq + n1 / sigma1_sq)
  mu1_post <- sigma1_post * (mu0 / tau0_sq + sum(y1) / sigma1_sq)
  theta1 <- rnorm(1, mean = mu1_post, sd = sqrt(sigma1_post))
  
  sigma2_post <- 1 / (1 / tau0_sq + n2 / sigma2_sq)
  mu2_post <- sigma2_post * (mu0 / tau0_sq + sum(y2) / sigma2_sq)
  theta2 <- rnorm(1, mean = mu2_post, sd = sqrt(sigma2_post))
  
  # Step 4: Update sigma1^2 and sigma2^2
  nu1_post <- nu0 + n1
  scale1_post <- (nu0 * sigma0_sq + sum((y1 - theta1)^2)) / nu1_post
  sigma1_sq <- rinvgamma(1, shape = nu1_post / 2, scale = nu1_post * scale1_post / 2)
  
  nu2_post <- nu0 + n2
  scale2_post <- (nu0 * sigma0_sq + sum((y2 - theta2)^2)) / nu2_post
  sigma2_sq <- rinvgamma(1, shape = nu2_post / 2, scale = nu2_post * scale2_post / 2)
  
  # Store samples
  theta1_samples[s] <- theta1
  theta2_samples[s] <- theta2
  p_samples[s] <- p
  sigma1_sq_samples[s] <- sigma1_sq
  sigma2_sq_samples[s] <- sigma2_sq
  
  # Step 5: Calculate theta_(1) and theta_(2) for this iteration
  theta_min_samples[s] <- min(theta1, theta2)
  theta_max_samples[s] <- max(theta1, theta2)
}

# Autocorrelation and effective sample sizes
acf(theta_min_samples, main = "Autocorrelation of Theta_min")
acf(theta_max_samples, main = "Autocorrelation of Theta_max")

effective_size_min <- effectiveSize(theta_min_samples)
effective_size_max <- effectiveSize(theta_max_samples)

cat("Effective Sample Size for Theta_min:", effective_size_min, "\n")
cat("Effective Sample Size for Theta_max:", effective_size_max, "\n")

# Compute and plot the autocorrelation functions for theta_min and theta_max
acf(theta_min_samples, main = "Autocorrelation of theta_min_samples")
acf(theta_max_samples, main = "Autocorrelation of theta_max_samples")

# Calculate effective sample sizes
theta_min_ess <- effectiveSize(theta_min_samples)
theta_max_ess <- effectiveSize(theta_max_samples)

print(paste("Effective Sample Size for theta_min:", theta_min_ess))
print(paste("Effective Sample Size for theta_max:", theta_max_ess))



"d) For each iteration s of the Gibbs sampler, sample a value x ∼
binary(p^(s)), then sample Y˜(s) ∼ normal(θx^(s), σx^2(s). 
Plot a histogram or kernel density estimate for the empirical distribution of
Y˜^(1), . . . , Y˜^(S), and compare to the distribution in part a). 
Discuss the adequacy of this two-component mixture model for the glucose
data."

# Number of iterations
num_iter <- length(p_samples)

# Store the simulated Y_tilde values
Y_tilde_samples <- numeric(num_iter)

# Simulate Y_tilde for each iteration
for (s in 1:num_iter) {
  # Sample group indicator x based on p_samples[s]
  x <- rbinom(1, 1, p_samples[s]) + 1
  
  # Sample Y_tilde based on the group parameters for this iteration
  if (x == 1) {
    Y_tilde_samples[s] <- rnorm(1, mean = theta1_samples[s], sd = sqrt(sigma1_sq_samples[s]))
  } else {
    Y_tilde_samples[s] <- rnorm(1, mean = theta2_samples[s], sd = sqrt(sigma2_sq_samples[s]))
  }
}

# Plot the histogram of Y_tilde_samples
hist(Y_tilde_samples, breaks = 30, probability = TRUE, 
     main = "Histogram of Simulated Y_tilde Values", 
     xlab = "Simulated Glucose Concentrations (Y_tilde)", 
     col = "lightblue")

# Overlay a kernel density estimate for Y_tilde
lines(density(Y_tilde_samples), col = "darkblue", lwd = 2)

# Overlay the original glucose data distribution from part (a) for comparison
lines(density(glucose_data$V1), col = "red", lwd = 2, lty = 2)

legend("topright", legend = c("Simulated Y_tilde", "Original Glucose Data"),
       col = c("darkblue", "red"), lwd = 2, lty = c(1, 2))