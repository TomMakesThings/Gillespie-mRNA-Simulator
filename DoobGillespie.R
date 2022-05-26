# Code to run the Doob-Gillespie algorithm

library(SciViews)
library(ggplot2)
library(stringr)
library(ggpubr)
library(scales)
library(reshape)
library(modi)

# Set seed for reproducability
set.seed(123)
# Sed directory
setwd("C:/Users/redds/Documents/GitHub/Gillespie")

func <- "B"

# Set parameters for function A
if (toupper(func) == "A") {
  # Set rates for production, degradation and complex formation
  rate_parameters <- list(lambda_1 = 20, beta_1 = 1, lambda_2 = 1, beta_2 = 0.1)
  
  # Rate formulas per reaction
  # Values are calculated using parameters specified in the gillespie function
  reaction_rates = c("rate_constants$lambda_1", # mRNA (x1) transcription
                     "rate_constants$beta_1 * X[[\"x1\"]]", # mRNA (x1) degradation
                     "rate_constants$lambda_2 * X[[\"x1\"]]", # x2 translation
                     "rate_constants$beta_2 * X[[\"x2\"]]", # x2 degradation
                     "rate_constants$lambda_2 * X[[\"x1\"]]", # x3 translation
                     "rate_constants$beta_2 * X[[\"x3\"]]") # x3 degradation
  
} else {
  # Set parameters for function B
  rate_parameters <- list(lambda_1 = 30, K = 40, beta_1 = 1, lambda_2 = 1, beta_2 = 0.1)
  
  reaction_rates = c("rate_constants$lambda_1 * (rate_constants$K) / (rate_constants$K + X[[\"x1\"]])", # mRNA (x1) transcription
                     "rate_constants$beta_1 * X[[\"x1\"]]", # mRNA (x1) degradation
                     "rate_constants$lambda_2 * X[[\"x1\"]]", # x2 translation
                     "rate_constants$beta_2 * X[[\"x2\"]]", # x2 degradation
                     "rate_constants$lambda_2 * X[[\"x1\"]]", # x3 translation
                     "rate_constants$beta_2 * X[[\"x3\"]]") # x3 degradation
}

# Set initial number of molecules
X_initial <- list(x1 = 50, x2 = 50, x3 = 50)

# Change in molecule numbers per reaction
reaction_deltas <- list(c(x1 = 1, x2 = 0, x3 = 0), c(x1 = -1, x2 = 0, x3 = 0),
                        c(x1 = 0, x2 = 1, x3 = 0), c(x1 = 0, x2 = -1, x3 = 0),
                        c(x1 = 0, x2 = 0, x3 = 1), c(x1 = 0, x2 = 0, x3 = -1))

# Dynamically calculate the reaction rate given an equation
calculateRate <- function(rate_equation, rate_constants, X) {
  # Converts string to equation and evaluates it
  # E.g. if variables beta = 10 and x1 = 5, then "beta * X$x1" becomes beta * x1 and 50 is returned
  eval(parse(text = rate_equation))
}

# Generate random wait time between reactions
calculateJumpTime <- function(rT) {
  # Generate random number between 0 - 1 from a uniform distribution
  u_time <- runif(1, 0, 1)
  # Calculate the random jump time by converting from a uniform to exponentially distributed number
  jump_time <- - ln(u_time) / rT
  
  return(jump_time)
}

# Calculate normalised variation, covariance and whether stationarity is reached
calculateStats <- function(X_over_time, rate_constants, jump_times) {
  CV_norm <- NA
  covar_norm <- NA
  
  # Extract the rate constant parameters
  lambda_1 <- rate_constants$lambda_1
  K <- rate_constants$K
  beta_1 <- rate_constants$beta_1
  lambda_2 <- rate_constants$lambda_2
  beta_2 <- rate_constants$beta_2
  
  # Calculate the average number of x1, x2 and x3 molecules
  x_averages <- list()
  
  for (x in rownames(X_over_time)) {
    x_averages[[x]] <- weighted.mean(x = unlist(X_over_time[x,]), w = jump_times)
  }
  
  # Calculate normalised variance
  eta_11 <- 1 / x_averages[["x1"]]
  
  # Normalised co-variance
  eta_12 <- eta_11 + (beta_2 / (beta_1 + beta_2))
  
  # Record results
  results <- list(eta_11 = eta_11,
                  eta_12 = eta_12)
  
  return(results)
}

# Doob-Gillespie algorithm
# X: initial number of molecules
# rate_constants: list of constant parameters for calculating rates
# deltas: list of vectors for the number of molecules that will change per reaction
# rates: formulas for reaction rates per reaction
# N: max number of iterations
# stop_at_stationary: whether to stop the simulation early if reached stationarity condition
# min_N: minimum number of iterations is stopping when reached stationarity
# flux_threshold: minimum difference between <R+> and <R-> for stationarity
# verbose: set as > 0 for debugging messages
gillespie <- function(X, rate_constants, deltas, rates, N = 100000,
                      stop_at_stationary = FALSE, min_N = 10000, 
                      flux_threshold = 0.1, verbose = 0) {
  # Record number of molecules over time in a matrix
  X_over_time <- matrix(data = X, nrow = length(X))
  # Record which reaction occurred at each time step t
  reaction_trace <- c()
  # Record time t between reactions
  times <- c(0)
  
  # Determine which reactions contains molecule degradation
  degradation_reactions <- which(unlist(lapply(deltas, function(x) any(x < 0))))
  
  # Calculate the magnitude of the highest degradation delta value
  degradation_deltas <- unlist(deltas[degradation_reactions])
  max_degration <- abs(degradation_deltas[which.max(abs(degradation_deltas))])
  
  n <- 0
  stationarity_reached <- FALSE
  
  # Run simulation N times or until stationarity is reached
  while (n < N & !stationarity_reached) {
    n <- n + 1
    
    # Get the number of molecules at time step t
    current_molecules <- unlist(X_over_time[, dim(X_over_time)[2]])
    names(current_molecules) <- names(X)
    
    # Initially set all reactions as possible at the current time step
    possible_reactions <- 1:length(rates)
    
    # Check if any degradation reactions will cause a negative number of molecules
    if (any(current_molecules - rep(max_degration, length(current_molecules)) < 0)) {
      # Find which reactions cannot occur
      impossible_reactions <- c()
      
      for (reaction_idx in degradation_reactions) {
        # Calculate number of molecules left if reaction occurs
        n_molecules <- unlist(X_over_time[, dim(X_over_time)[2]]) + deltas[[reaction_idx]]
        
        # Record reaction as impossible if the number of any molecule is below zero
        if (any(n_molecules < 0)) {
          impossible_reactions <- c(impossible_reactions, reaction_idx)
        }
      }
      # Update list of possible reactions to remove impossible ones
      possible_reactions <- possible_reactions[!possible_reactions %in% impossible_reactions]
    }
    
    # Calculate the current reaction rates for the system at time t, excluding impossible reactions
    current_rates <- unname(unlist(lapply(rates[possible_reactions],
                                          function(rate_equation) calculateRate(rate_equation,
                                                                                rate_constants,
                                                                                current_molecules))))
    # Calculate the total reaction rate
    r_total <- sum(current_rates)
    
    # Generate the jump time for the next reaction
    times <- c(times, calculateJumpTime(r_total))
    
    # Calculate the probability of each reaction occurring
    reaction_probs <- current_rates / r_total
    prob_partitions <- cumsum(reaction_probs)

    # Select a random reaction (relative to the probability of each occurring)
    u_reaction <- runif(1, 0, 1)
    reaction_idx <- min(which(prob_partitions > u_reaction))
    
    # Calculate the molecule numbers after the chosen reaction
    n_molecules <- current_molecules + deltas[possible_reactions][[reaction_idx]]
    
    # Record which reaction occurred
    reaction_trace <- c(reaction_trace, reaction_idx)
    # Update the number of molecules by the chosen reaction
    X_over_time <- cbind(X_over_time, n_molecules)
    
    if (verbose > 0 & n %% 10000 == 0) {
      message(paste("Reached iteration", n))
    }
    
    if (stop_at_stationary & n %% 500 == 0 & n >= min_N) {
      # Check if difference in birth and death flux is less than a threshold
      stats <- calculateStats(X_over_time, rate_constants, times,
                              flux_threshold = flux_threshold,
                              calculate_efficiency = FALSE,
                              calculate_variation = FALSE,
                              calculate_covariance = FALSE)
      
      # Stop simulation if TRUE
      stationarity_reached <- stats$stationarity_reached
      # flux_differences <- stats$flux_differences
      
      if (stationarity_reached) {
        if (verbose > 0) {
          message(paste("Stationarity reached after", n, "iterations"))
        }
      }
    }
  }
  
  # Rename columns to time t for each reaction
  colnames(X_over_time) <- cumsum(times)
  
  return(list(molecules_over_time = X_over_time,
              reaction_trace = reaction_trace,
              jump_times = times,
              n_iterations = n))
}

# Run the Gillespie algorithm
simulation_results <- gillespie(X = X_initial,
                                rate_constants = rate_parameters,
                                deltas = reaction_deltas,
                                rates = reaction_rates,
                                N = 250000,
                                stop_at_stationary = FALSE,
                                min_N = 1000,
                                verbose = 1)

message(signif(mean(unlist(simulation_results$molecules_over_time["x1",])), 5))

# saveRDS(simulation_results, file = "sim_A.rds")
# saveRDS(simulation_results, file = "sim_B.rds")

# Open the saved simulations for functions (a) and (b)
sim_A <- readRDS(file = "sim_A.rds")
sim_B <- readRDS(file = "sim_B.rds")

# Calculate efficiency for the simulation
simulation_stats <- calculateStats(X_over_time = simulation_results$molecules_over_time,
                                   rate_constants = rate_parameters,
                                   jump_times = simulation_results$jump_times)

# Create plot of molecule abundance over time
timePlot <- function(molecule_counts, plot_colours = c("red", "deepskyblue", "purple"),
                     title = "Number of Molecules Over Time",
                     xaxis_label = "Time (s)", yaxis_label = "Number of Molecules",
                     xlimits = NA, ylimits = NA) {
  # Convert data into form for plotting
  molecules_time_df <- data.frame(n_molecules = unlist(c(molecule_counts)),
                                  molecule_type = row.names(molecule_counts),
                                  time = as.numeric(rep(colnames(molecule_counts),
                                                        each = nrow(molecule_counts))))
  # Plot molecules over time
  time_plot <- ggplot(data = molecules_time_df,
                      aes(x = time, y = n_molecules,
                          group = molecule_type,
                          color = molecule_type)) +
    geom_line() +
    labs(title = title,
         x = xaxis_label,
         y = yaxis_label) +
    guides(color = guide_legend(title = "Molecule Type")) +
    scale_color_manual(labels = c(expression("mRNA (" ~ x[1] ~ ")"),
                                  expression("Protein (" ~ x[2] ~ ")"),
                                  expression("Protein (" ~ x[3] ~ ")")),
                       values = plot_colours)
  
  if (!is.na(xlimits)) {
    # Set x-axis limits
    time_plot <- time_plot + scale_x_continuous(limits = xlimits)
  }
  if (!is.na(ylimits)) {
    # Set y-axis limits
    time_plot <- time_plot + scale_y_continuous(limits = ylimits)
  }
  
  return(time_plot)
}

# Create plots of molecules over time
time_plot_A <- timePlot(sim_A$molecules_over_time, 
                        title = "Function A",
                        xlimits = c(0, 2100),
                        ylimits = c(0, 300))

time_plot_B <- timePlot(sim_B$molecules_over_time,
                        title = "Function B",
                        xlimits = c(0, 2100),
                        ylimits = c(0, 300))

# Save as subplots
pdf("Molecules_Over_Time.pdf", width = 18, height = 8)
time_subplots <- ggarrange(time_plot_A , time_plot_B, nrow = 2, widths = c(1, 1),
                           common.legend = TRUE, legend = "right")
annotate_figure(time_subplots, top = text_grob("Molecules Over Time", 
                                               face = "bold", size = 14))
dev.off()


# Test different parameters
parameterSearch <- function(test_lambdas_1 = c(1), test_lambdas_2 = c(1),
                            test_betas_1 = c(1), test_betas_2 = c(1),
                            test_Ks = c(NA), X, deltas, rates, N = 50000,
                            flux_threshold = 0.1, stop_at_stationary = FALSE,
                            min_N = 10000, verbose = 0) {
  
  # Record parameters and their results
  parameter_combination <- data.frame()
  sim_results <- list()
  stats_results <- list()
  
  # Assign an ID for the combination
  id <- 0
  
  # Test each parameter combination
  for (l1 in test_lambdas_1) {
    for (l2 in test_lambdas_2) {
      for (b1 in test_betas_1) {
        for (b2 in test_betas_2) {
          for (k in test_Ks) {
            # Set unique ID
            id <- id + 1
            
            new_params <- list(lambda_1 = l1,
                               lambda_2 = l2,
                               beta_1 = b1,
                               beta_2 = b2,
                               K = k)
            
            # Record the lambdas, betas and K
            parameter_combination <- rbind(parameter_combination, new_params)
            
            # Run the Gillespie algorithm
            test_simulation <- gillespie(X = X_initial,
                                         rate_constants = new_params,
                                         deltas = deltas,
                                         rates = rates,
                                         N = N,
                                         flux_threshold = flux_threshold,
                                         stop_at_stationary = stop_at_stationary,
                                         min_N = min_N,
                                         verbose = verbose)
            
            # Calculate normalised variance and covariances
            stats <- calculateStats(X_over_time = test_simulation$molecules_over_time,
                                    rate_constants = list(lambda = l, beta = b, C = c),
                                    jump_times = test_simulation$jump_times)
            
            # Record the results for the parameter combination
            sim_results[[id]] <- test_simulation
            stats_results[[id]] <- stats
          }
        }
      }
    }
  }
  
  return(list(sim_results = sim_results,
              stats_results = stats_results,
              parameter_combination = parameter_combination))
}

# Parameters to test
test_lambdas_1 <- c(10, 20, 30, 40, 50, 60, 70, 80)
test_Ks <- test_lambdas_1

# Test different combinations of parameters and return the results of each
search_results <- parameterSearch(test_lambdas = test_lambdas, 
                                  test_betas = test_betas,
                                  test_Cs = test_Cs,
                                  X = X_initial,
                                  deltas = reaction_deltas,
                                  rates = reaction_rates,
                                  N = 1000)

setwd("C:/Users/redds/Documents/GitHub/SystemsBiologyGroup")

# Save results to file
# saveRDS(search_results, file = "parameter_test_results.rds")
search_results <- readRDS("parameter_test_results_100000_final.rds")

# Get the simulation time trace
parameter_sims <- search_results$sim_results
# Combine lambda, beta, C data with stats
parameter_stats_df <- cbind(search_results$parameter_combination,
                            do.call(rbind, search_results$stats_results))
# Convert columns from list to numeric and logical
parameter_stats_df[c(4:7)] <- lapply(parameter_stats_df[c(4:7)], as.numeric)
parameter_stats_df[c(8)] <- lapply(parameter_stats_df[c(8)], as.logical)
# Covert C to an ordered factor
parameter_stats_df$C <- factor(parameter_stats_df$C, levels = unique(parameter_stats_df$C))

# Calculate normalised standard deviation
parameter_stats_df$normalised_sd <- sqrt(parameter_stats_df$normalised_variance)

# Calculate the absolute mean difference between fluxes <R+> - <R-> for each molecule
parameter_stats_df$mean_flux_diff <- abs(unlist(lapply(parameter_stats_df$flux_differences,
                                                       function(x) mean(unlist(x)))))