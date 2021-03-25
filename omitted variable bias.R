#------------------------ libraries ------------------------
library(tidyverse)
library(tidylog)
library(MASS)
library(ggthemes)

set.seed(5678)

#------------------------ parameters ------------------------
betas <- matrix(c(-3, 2, 0.5), 
                nrow = 3, 
                ncol = 1)

alphas <-  matrix(c(1.5, -3), 
                 nrow = 2, 
                 ncol = 1)

sigmas <- matrix(c(1, 0.75), 
                 nrow = 2, 
                 ncol = 1)

#------------------------ functions ------------------------
pseudo_data_generator <- function(n_obs, betas, alphas, sigmas) {
  X_i <- as.matrix(rnorm(1:n_obs, mean = 10, sd = 5))
  
  error_1 <- as.matrix(rnorm(1:n_obs), mean = 0, sd = sigmas[1, 1])
  error_2 <- as.matrix(rnorm(1:n_obs), mean = 0, sd = sigmas[2, 1])
  
  D <- alphas[1,1] + alphas[2,1]*X_i + error_2
  
  y <- betas[1,1] + betas[2,1]*X_i + betas[3,1]*D + error_1
  
  return(list("X_i" = X_i, 
              "D" = D, 
              "Y" = y))
}

myOLS <- function(data, outcomes) {
  Y <- outcomes
  X <- data
  
  degrees_of_freedom <- nrow(X)-ncol(X)
  
  beta_estimates <- ginv(t(X)%*%X)%*%t(X)%*%Y
  
  error_estimates   <- Y - X%*%beta_estimates
  sigma_estimate    <- (t(error_estimates)%*%error_estimates) / degrees_of_freedom
  variance_estimate <- ginv(t(X)%*%X)
  standard_error    <- as.matrix(sqrt(diag(variance_estimate)))
  
  RSS       <- t(error_estimates)%*%error_estimates         
  y_bar     <- mean(Y)
  yybar     <- Y-y_bar
  TSS       <- t(yybar)%*%yybar 
  R_squared <- 1 - RSS/TSS 
  
  t_statistic <- beta_estimates / standard_error
  p_values     <- 2*(1-pt(abs(t_statistic), degrees_of_freedom)) %>%
                 round(5)
  
  estimates <- tibble("estimate"    = c(beta_estimates), 
                      "std_error"   = c(standard_error), 
                      "t_statistic" = c(t_statistic), 
                      "p_value"     = c(p_values))
  
  list("estimates" = estimates,
       "r_squared" = R_squared)
}

omit_variable <- function(n_obs, betas, alphas, sigmas, omit = TRUE) {
  pseudo_data <- pseudo_data_generator(n_obs, betas, alphas, sigmas)
  outcomes <- pseudo_data$Y
  
  if(omit == FALSE) {
    data <- cbind(rep(1, n_obs), pseudo_data$X_i, pseudo_data$D)
  } 
  
  if(omit == TRUE) {
    data <- cbind(rep(1, n_obs), pseudo_data$X_i)
  }
  
  results <- myOLS(data, outcomes)[[1]]
  
  return(as.matrix(results))
}

clean_results <- function(output, omit = TRUE) {
  results <- c()
  
  for(i in 1:dim(output)[3]) {
    results <- rbind(results, output[,,i])
  }
  
  if(omit == FALSE) {
    beta_0 <- results[seq(1,(dim(output)[3]*2), 3),]
    beta_1 <- results[seq(2,(dim(output)[3]*2), 3),]
  }
  
  if(omit == TRUE) {
    beta_0 <- results[seq(1,(dim(output)[3]*2), 2),]
    beta_1 <- results[seq(2,(dim(output)[3]*2), 2),]
  }
  
  return(list("beta_0" = as_tibble(beta_0), 
              "beta_1" = as_tibble(beta_1)))
}

#------------------------ analysis ------------------------
n_obs <- 100

pseudo_data <- pseudo_data_generator(n_obs, betas, alphas, sigmas)
data <- cbind(rep(1, n_obs), pseudo_data$X_i, pseudo_data$D)
outcomes <- pseudo_data$Y

myOLS(data, outcomes)[[1]] %>% names()

#------------------------ Monte Carlo ------------------------
n_obs <- 1000

output <- replicate(1000, omit_variable(n_obs, betas, alphas, sigmas, omit = T))

results <- clean_results(output, omit = T)

results <- bind_rows(results, .id = "parameter")

results1 <- array(numeric(), c(3, 4, 100))

for(i in 1:100) {
  pseudo_data <- pseudo_data_generator(n_obs, betas, alphas, sigmas)
  outcomes <- pseudo_data$Y
  
  if(omit == FALSE) {
    data <- cbind(rep(1, n_obs), pseudo_data$X_i, pseudo_data$D)
  } 
  
  if(omit == TRUE) {
    data <- cbind(rep(1, n_obs), pseudo_data$X_i)
  }
  
  results <- myOLS(data, outcomes)[[1]]
  
  
}


results %>%
  ggplot() +
  geom_histogram(aes(estimate), 
                 fill = "royalblue", 
                 color = "navyblue", 
                 bins = 30) +
  #geom_vline(xintercept = betas[1], color = "red") +
  facet_wrap(.~parameter, scales = "free") +
  theme_linedraw() +
  theme(axis.title.y = element_blank())

a#------------------------ Bias ------------------------
n_obs <- 1000

# theoretical bias
pseudo_data <- pseudo_data_generator(n_obs, betas, alphas, sigmas)
estimate <- myOLS(pseudo_data$X_i, pseudo_data$D)[[1]][1]

theoretical_bias <- estimate*betas[3]


