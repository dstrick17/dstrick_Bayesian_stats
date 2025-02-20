"7.6 Diabetes data: A population of 532 women living near Phoenix, Arizona were 
tested for diabetes. Other information was gathered from these
women at the time of testing, including number of pregnancies, glucose
level, blood pressure, skin fold thickness, body mass index, diabetes pedigree 
and age. This information appears in the file azdiabetes.dat. Model
the joint distribution of these variables for the diabetics and non-diabetics
separately, using a multivariate normal distribution:"


# Load the glucose data from GitHub
azdiabetes_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrickMA578/refs/heads/main/CASMA578/HW4/azdiabetes.txt", header=TRUE)
# Check the first few rows to confirm the structure of the data
head(azdiabetes_data)
library(MCMCpack)

"a) For both groups separately, use the following type of unit information
prior, where Σˆ is the sample covariance matrix.
i. µ0 = y¯, Λ0 = Σˆ;
ii. S0 = Σˆ, ν0 = p + 2 = 9 .
Generate at least 10,000 Monte Carlo samples for {θd, Σd} and
{θn, Σn}, the model parameters for diabetics and non-diabetics respectively. For each of the seven variables j ∈ {1, . . . , 7}, compare the
marginal posterior distributions of θd,j and θn,j . Which variables seem
to differ between the two groups? Also obtain Pr(θd,j > θn,j |Y) for
each j ∈ {1, . . . , 7}."

library(dplyr)

# Split data into diabetic and non-diabetic groups
library(dplyr)
diabetic_group <- azdiabetes_data %>% filter(diabetes == "Yes")
non_diabetic_group <- azdiabetes_data %>% filter(diabetes == "No")

head(diabetic_group)
head(non_diabetic_group)

# Calculate means and covariances for each group
mu_d <- colMeans(diabetic_group[, 1:7])  # Assuming first 7 columns are variables
mu_n <- colMeans(non_diabetic_group[, 1:7])

cov_d <- cov(diabetic_group[, 1:7])
cov_n <- cov(non_diabetic_group[, 1:7])

# Set priors based on sample mean and covariance
nu0 <- 9  # p + 2 = 7 + 2
S0_d <- cov_d
S0_n <- cov_n

# Gibbs sampler parameters
num_iter <- 10000
theta_d_samples <- matrix(0, nrow = num_iter, ncol = 7)
theta_n_samples <- matrix(0, nrow = num_iter, ncol = 7)

# Add storage for Sigma_d and Sigma_n samples
sigma_d_samples <- array(0, dim = c(7, 7, num_iter))  # 7x7 matrix for each iteration
sigma_n_samples <- array(0, dim = c(7, 7, num_iter))  # 7x7 matrix for each iteration



# Gibbs sampling loop
for (i in 1:num_iter) {
  # Sample covariance matrix from inverse-Wishart
  Sigma_d <- riwish(nu0, S0_d)
  Sigma_n <- riwish(nu0, S0_n)
  
  # Sample mean vector from multivariate normal
  theta_d <- mvrnorm(1, mu = mu_d, Sigma = Sigma_d / nrow(diabetic_group))
  theta_n <- mvrnorm(1, mu = mu_n, Sigma = Sigma_n / nrow(non_diabetic_group))
  
  # Store samples
  theta_d_samples[i, ] <- theta_d
  theta_n_samples[i, ] <- theta_n
  sigma_d_samples[, , i] <- Sigma_d
  sigma_n_samples[, , i] <- Sigma_n
}

# Calculate Pr(θd,j > θn,j | Y) for each variable j
prob_greater <- colMeans(theta_d_samples > theta_n_samples)
print(prob_greater)

###############################################################################3
"b) Obtain the posterior means of Σd and Σn, and plot the entries versus
each other. What are the main differences, if any?"

# Variable names
variable_names <- c("npreg", "glu", "bp", "skin", "bmi", "ped", "age")

# Calculate posterior means for Sigma_d and Sigma_n
posterior_mean_Sigma_d <- apply(sigma_d_samples, c(1, 2), mean)
posterior_mean_Sigma_n <- apply(sigma_n_samples, c(1, 2), mean)

# Assign row and column names to the matrices
rownames(posterior_mean_Sigma_d) <- variable_names
colnames(posterior_mean_Sigma_d) <- variable_names
rownames(posterior_mean_Sigma_n) <- variable_names
colnames(posterior_mean_Sigma_n) <- variable_names

# Print posterior means for Σd and Σn
cat("Posterior mean of Σd (Diabetic group):\n")
print(posterior_mean_Sigma_d)
cat("\nPosterior mean of Σn (Non-diabetic group):\n")
print(posterior_mean_Sigma_n)

# Flatten the matrices for comparison
sigma_d_entries <- as.vector(posterior_mean_Sigma_d)
sigma_n_entries <- as.vector(posterior_mean_Sigma_n)

# Create labels based on variable pairs
labels <- outer(variable_names, variable_names, paste, sep = "_")
labels <- as.vector(labels)

# Create a data frame for plotting with labels
comparison_df <- data.frame(Sigma_d = sigma_d_entries, Sigma_n = sigma_n_entries, Label = labels)

# Plot with ggplot2 and add labels for each point
library(ggplot2)
ggplot(comparison_df, aes(x = Sigma_d, y = Sigma_n)) +
  geom_point() +
  geom_text(aes(label = Label), hjust = 1.2, vjust = 1.2, size = 3) +  # Add labels
  labs(
    title = "Posterior Mean Entries of Covariance Matrices Σd and Σn",
    x = expression(Sigma[d]),
    y = expression(Sigma[n])
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal()