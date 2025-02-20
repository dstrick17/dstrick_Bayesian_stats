# 9.1 Extrapolation: The file swim.dat contains data on the amount of time,
# in seconds, it takes each of four high school swimmers to swim 50 yards.
# Each swimmer has six times, taken on a biweekly basis.

# a) Perform the following data analysis for each swimmer separately:
    # i. Fit a linear regression model of swimming time as the response and
    # week as the explanatory variable. To formulate your prior, use the
    # information that competitive times for this age group generally
    # range from 22 to 24 seconds

    # ii. For each swimmer j, obtain a posterior predictive distribution
    # for Y∗j, their time if they were to swim two weeks from the last
    # recorded time.

# b) The coach of the team has to decide which of the four swimmers will
# compete in a swimming meet in two weeks. Using your predictive distributions, 
# compute Pr(Y∗j = max{Y ∗ 1, . . . , Y ∗4 }|Y)) for each swimmer j, and based
# on this make a recommendation to the coach.
# Compute probabilities of winning

# Load data set
swim_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrick_Bayesian_stats/refs/heads/main/CASMA578/HW5/swim.txt", header = FALSE)
#
# View the data structure
print(swim_data)

library(dplyr)
library(rstanarm)
library(bayestestR)
# Add a Swimmer column to identify each swimmer
swim_data$Swimmer <- 1:nrow(swim_data)

# i. Fit Bayesian linear regression models for each swimmer
fit_swimmer_model <- function(swimmer_times) {
  # Prepare data for the model
  swimmer_data <- data.frame(
    Week = 1:6,
    Time = as.numeric(swimmer_times)
  )
  
  # Fit Bayesian linear regression model
  model <- stan_glm(
    Time ~ Week,
    data = swimmer_data,
    prior_intercept = normal(23, 1),  # Prior for intercept based on 22-24 second range
    prior = normal(0, 0.1),           # Prior for slope (expecting minimal change)
    prior_aux = exponential(1),       # Prior for residual standard deviation
    seed = 123                        # For reproducibility
  )
  
  return(model)
}

# Part (a) Detailed Analysis
cat("PART (a): Detailed Swimmer Analysis\n")
cat("==================================\n\n")

# Stores results for detailed output
swimmers_analysis <- list()

# Fit models for each swimmer
bayesian_results <- list()
for (i in 1:nrow(swim_data)) {
  # Fit model for this swimmer
  swimmer_model <- fit_swimmer_model(swim_data[i, 1:6])
  
  # Store the model
  bayesian_results[[i]] <- swimmer_model
  
  # Prepare detailed analysis for this swimmer
  swimmers_analysis[[i]] <- list(
    original_times = as.numeric(swim_data[i, 1:6]),
    model_summary = summary(swimmer_model)
  )
  
  # i. Detailed Model Output
  cat("i. Bayesian Regression Results for Swimmer", i, ":\n")
  print(summary(swimmer_model))
  cat("\n")
}

# ii. Posterior Predictive Distribution
cat("ii. Posterior Predictive Distributions\n")
cat("-------------------------------------\n")

# Predict times for each swimmer two weeks from last recorded time
prediction_results <- list()
for (i in 1:length(bayesian_results)) {
  # Prepare prediction data
  prediction_data <- data.frame(Week = 7)
  
  # Generate posterior predictions
  predictions <- posterior_predict(
    bayesian_results[[i]], 
    newdata = prediction_data,
    draws = 1000
  )
  
  # Store predictions
  prediction_results[[i]] <- predictions
  
  # Print detailed distribution for each swimmer
  cat("Swimmer", i, "Predicted Time Distribution:\n")
  print(distribution_normal(predictions))
  cat("Predicted Time (Mean ± SD):", 
      round(mean(predictions), 2), 
      "±", 
      round(sd(predictions), 2), 
      "seconds\n\n")
}

# Part (b) Winning Probabilities
cat("PART (b): Competitive Probability Analysis\n")
cat("=========================================\n")

# Compute probabilities of winning
win_probabilities <- numeric(length(bayesian_results))
for (j in 1:length(bayesian_results)) {
  # Compute probability of this swimmer being the fastest
  win_probabilities[j] <- mean(
    apply(sapply(prediction_results, function(x) x[,1]), 1, 
          function(swimmer_times) all(swimmer_times[j] == min(swimmer_times)))
  )
}

# Print detailed win probabilities
cat("Win Probabilities for Each Swimmer:\n")
for (j in 1:length(win_probabilities)) {
  cat("Swimmer", j, ":", 
      round(win_probabilities[j] * 100, 2), 
      "% chance of being the fastest\n")
}

# Recommend the swimmer with the highest win probability
best_swimmer <- which.max(win_probabilities)
cat("\nCoach's Recommendation:\n")
cat("Swimmer", best_swimmer, 
    "has the highest chance of winning with a probability of", 
    round(max(win_probabilities) * 100, 2), 
    "%\n")

print(swim_data)

win_prob_df <- data.frame(
  Swimmer = paste("Swimmer", 1:length(win_probabilities)),
  Probability = win_probabilities
)

library(ggplot2)

for (i in 1:nrow(swim_data)) {
  swimmer_data <- data.frame(
    Week = 1:6,
    Time = as.numeric(swim_data[i, 1:6])
  )
  
  # Extract fitted regression line
  fit_line <- data.frame(
    Week = 1:7,
    Time = posterior_epred(bayesian_results[[i]], newdata = data.frame(Week = 1:7)) %>% colMeans()
  )}
  
  
  ########################### Plots
  library(dplyr)
  
  # Ensure required libraries are loaded
  library(ggplot2)
  library(dplyr)
  library(rstanarm)
  
  # Loop for plotting swimmer data and regression fit
  for (i in 1:nrow(swim_data)) {
    # Create swimmer data
    swimmer_data <- data.frame(
      Week = 1:6,
      Time = as.numeric(swim_data[i, 1:6])
    )
    
    # Ensure model exists for the swimmer
    if (!is.null(bayesian_results[[i]])) {
      # Generate predictions for the fitted line
      fit_line <- data.frame(
        Week = 1:7,
        Time = colMeans(posterior_epred(bayesian_results[[i]], newdata = data.frame(Week = 1:7)))
      )
      
      # Plot swimmer's data and regression fit
      print(
        ggplot(swimmer_data, aes(x = Week, y = Time)) +
          geom_point(color = "blue", size = 3) +
          geom_line(data = fit_line, aes(x = Week, y = Time), color = "red") +
          labs(title = paste("Swimmer", i, "- Bayesian Regression Fit"),
               x = "Week",
               y = "Time (seconds)") +
          theme_minimal()
      )
    } else {
      cat("Model for swimmer", i, "not found or not computed.\n")
    }
  }
  
  # Loop for plotting predicted time distributions
  for (i in 1:length(prediction_results)) {
    print(
      ggplot(data.frame(PredictedTime = prediction_results[[i]][, 1]), aes(x = PredictedTime)) +
        geom_density(fill = "lightblue", alpha = 0.6) +
        labs(title = paste("Swimmer", i, "- Predicted Time Distribution"),
             x = "Predicted Time (seconds)",
             y = "Density") +
        theme_minimal()
    )
  }
  
  # Prepare win probabilities dataframe
  win_prob_df <- data.frame(
    Swimmer = paste("Swimmer", 1:length(win_probabilities)),
    Probability = win_probabilities
  )
  
  # Plot win probabilities
  print(
    ggplot(win_prob_df, aes(x = Swimmer, y = Probability)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(title = "Win Probabilities for Each Swimmer",
           x = "Swimmer",
           y = "Probability of Being the Fastest") +
      theme_minimal()
  )
  