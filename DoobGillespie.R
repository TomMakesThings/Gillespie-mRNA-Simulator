# Code to run the Doob-Gillespie algorithm

library(SciViews)
library(ggplot2)
library(stringr)
library(ggpubr)

# Set seed for reproducability
set.seed(123)
# Set directory
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

# Calculate normalised variation and covariance
calculateStats <- function(X_over_time, rate_constants, jump_times,
                           function_type = "A") {
  
  # Extract the rate constant parameters
  lambda_1 <- rate_constants$lambda_1
  lambda_2 <- rate_constants$lambda_2
  beta_1 <- rate_constants$beta_1
  beta_2 <- rate_constants$beta_2
  K <- rate_constants$K
  
  # Calculate the average number of x1, x2 and x3 molecules
  x_averages <- list()
  
  for (x in rownames(X_over_time)) {
    x_averages[[x]] <- weighted.mean(x = unlist(X_over_time[x,]), w = jump_times)
  }
  
  # Calculate analytic etas
  if (toupper(function_type) == "A") {
    # Calculate normalised variance for function (a)
    eta_11 <- 1 / x_averages[["x1"]]
    
    # Normalised co-variance
    eta_12 <- eta_11 * (beta_2 / (beta_1 + beta_2))
    eta_23 <- eta_12
    
  } else {
    # Calculate placeholder alpha for function (b)
    alpha <- x_averages[["x1"]] / (K + x_averages[["x1"]])
    
    eta_11 <- 1 / (x_averages[["x1"]] + alpha)
    eta_12 <- eta_11 * (beta_2 / (beta_1 + (beta_1 * alpha) + beta_2))
    eta_23 <- eta_12
  }
  
  # Calculate weighted covariance between x1, x2 and x3 molecules
  abundance_df <- data.frame(t(X_over_time))
  abundance_df <- data.frame(lapply(abundance_df, as.numeric))
  molecule_cov <- cov.wt(abundance_df, as.vector(jump_times))$cov
  
  # Calculate numeric etas
  num_eta11 <- molecule_cov["x1", "x1"] / (x_averages[["x1"]] * x_averages[["x1"]])
  num_eta12 <- molecule_cov["x1", "x2"] / (x_averages[["x1"]] * x_averages[["x2"]])
  num_eta23 <- molecule_cov["x2", "x3"] / (x_averages[["x2"]] * x_averages[["x3"]])
  
  # Record results
  results <- list(analytical = list(eta_11 = eta_11,
                                    eta_12 = eta_12,
                                    eta_23 = eta_23),
                  numerical = list(eta_11 = num_eta11,
                                   eta_12 = num_eta12,
                                   eta_23 = num_eta23))
  
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

# saveRDS(simulation_results, file = "sim_A.rds")
# saveRDS(simulation_results, file = "sim_B.rds")

# Open the saved simulations for functions (a) and (b)
sim_A <- readRDS(file = "sim_A.rds")
sim_B <- readRDS(file = "sim_B.rds")

# Display average molecule abundance
message(paste("Average <x1>:", signif(mean(unlist(sim_A$molecules_over_time["x1",])), 5)))
message(paste("Average <x2>:", signif(mean(unlist(sim_A$molecules_over_time["x2",])), 5)))
message(paste("Average <x3>:", signif(mean(unlist(sim_A$molecules_over_time["x3",])), 5)))

message(paste("Average <x1>:", signif(mean(unlist(sim_B$molecules_over_time["x1",])), 5)))
message(paste("Average <x2>:", signif(mean(unlist(sim_B$molecules_over_time["x2",])), 5)))
message(paste("Average <x3>:", signif(mean(unlist(sim_B$molecules_over_time["x3",])), 5)))

# Calculate efficiency for the simulation
simulation_stats <- calculateStats(X_over_time = simulation_results$molecules_over_time,
                                   rate_constants = rate_parameters,
                                   jump_times = simulation_results$jump_times,
                                   function_type = func)

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
time_subplots <- ggarrange(time_plot_A, time_plot_B, nrow = 2, widths = c(1, 1),
                           common.legend = TRUE, legend = "right")
annotate_figure(time_subplots, top = text_grob("Molecules Over Time", 
                                               face = "bold", size = 14))
dev.off()

# Test different parameters
parameterSearch <- function(test_lambdas_1 = c(1), test_lambdas_2 = c(1),
                            test_betas_1 = c(1), test_betas_2 = c(1),
                            test_Ks = c(NA), X, deltas, rates, N = 50000,
                            flux_threshold = 0.1, stop_at_stationary = FALSE,
                            min_N = 10000, function_type = "A", verbose = 0,
                            previous_search_stats = data.frame()) {
  
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
            
            new_params <- list(lambda_1 = l1,
                               lambda_2 = l2,
                               beta_1 = b1,
                               beta_2 = b2,
                               K = k)
            
            # Check whether combination was tested previously in a different run
            continue_search <- TRUE
            
            if (nrow(previous_search_stats) > 0) {
              if (nrow(merge(new_params, previous_search_stats)) > 0) {
                continue_search <- FALSE
              }
            }
            
            if (continue_search) {
              # Set unique ID
              id <- id + 1
              
              if (verbose > 0) {
                message(paste("Reached step", id))
                message(paste("Testing", toString(new_params)))
              }
              
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
                                           min_N = min_N)
              
              # Calculate normalised variance and covariances
              stats <- calculateStats(X_over_time = test_simulation$molecules_over_time,
                                      rate_constants = new_params,
                                      jump_times = test_simulation$jump_times,
                                      function_type = function_type)
              
              # Record the results for the parameter combination
              sim_results[[id]] <- test_simulation
              stats_results[[id]] <- stats
            }
          }
        }
      }
    }
  }
  
  # Combine parameters with stats
  parameter_stats <- cbind(parameter_combination, do.call(rbind, stats_results))
  parameter_stats <- data.frame(lapply(parameter_stats, as.numeric))
  
  return(list(sim_results = sim_results,
              parameter_stats = parameter_stats))
}

if (toupper(func) == "A") {
  # Parameters to test
  test_lambdas_1 <- c((1:10)/10, (1:50)*2)
  test_Ks <- NA
} else {
  test_lambdas_1 <- c(0.1, 0.5, 1, 5, 10, 20, 50, 100)
  test_Ks <- test_lambdas_1
}

# Test different combinations of parameters and return the results of each
search_results <- parameterSearch(test_lambdas_1 = test_lambdas_1, 
                                  test_lambdas_2 = c(1),
                                  test_betas_1 = c(1),
                                  test_betas_2 = c(0.1),
                                  test_Ks = test_Ks,
                                  X = X_initial,
                                  deltas = reaction_deltas,
                                  rates = reaction_rates,
                                  N = 100000,
                                  function_type = "B",
                                  verbose = 1)

# saveRDS(search_B, file = "function_B_search2.rds")

# Open parameter search results for functions (a) and (b)
search_A <- readRDS(file = "function_A_search_new.rds")
search_B <- readRDS(file = "function_B_search_new.rds")

# ## Update incorrect stats for B
# params_B <- list()
# stats_results_B <- list()
# 
# for (r in 1:length(search_A$sim_results)) {
#   params_B[[r]] = list(lambda_1 = search_A$parameter_stats[r, "lambda_1"],
#                        lambda_2 = search_A$parameter_stats[r, "lambda_2"],
#                        beta_1 = search_A$parameter_stats[r, "beta_1"],
#                        beta_2 = search_A$parameter_stats[r, "beta_2"],
#                        K = search_A$parameter_stats[r, "K"])
#   stats_results_B[[r]] <- calculateStats(X_over_time = search_A$sim_results[[r]]$molecules_over_time,
#                                          rate_constants = params_B[[r]],
#                                          jump_times = search_A$sim_results[[r]]$jump_times,
#                                          function_type = "A")
# }
# 
# numeric_stats <- cbind(do.call(rbind.data.frame, params_B),
#                            do.call(rbind, lapply(stats_results_B, function(x) x$numerical)))
# numeric_stats <- data.frame(lapply(numeric_stats, as.numeric))
# analytic_stats <- cbind(do.call(rbind.data.frame, params_B),
#                         do.call(rbind, lapply(stats_results_B, function(x) x$analytical)))
# analytic_stats <- data.frame(lapply(analytic_stats, as.numeric))
# 
# search_A$numerical_stats <- numeric_stats
# search_A$analytical_stats <- analytic_stats
# 
# saveRDS(search_B, file = "function_B_search_new.rds")

# rps <- as.list(search_A$parameter_stats[1,][1:5])
# rps <- lapply(rps, as.numeric)
# 
# calculateStats(X_over_time = search_A$sim_results[[1]]$molecules_over_time,
#                rate_constants = rps,
#                jump_times = search_B$sim_results[[r]]$jump_times,
#                function_type = "A")


plotCovariance <- function(points_A, points_B, no_outliers_A, no_outliers_B) {
  # Combine eta values for each function into a dataframe
  combined_stats <- rbind(points_A, points_B)
  combined_stats$func <- c(rep("A", nrow(points_A)),
                           rep("B", nrow(points_B)))
  
  # Create dataframe for trend line in which outliers are removed
  no_outlier_stats <- rbind(no_outliers_A, no_outliers_B)
  no_outlier_stats$func <- c(rep("A", nrow(no_outliers_A)),
                             rep("B", nrow(no_outliers_B)))
  
  # Plot eta 11 and eta 12
  var_covar_plot <- ggplot(data = combined_stats, mapping = aes(x = eta_11, y = eta_12,
                                         color = func, shape = func)) +
    geom_point(data = combined_stats, size = 2.5, alpha = 0.5) +
    labs(title = expression("Normalized Variance " ~ eta[11] ~ "and Normalized Covariance" ~ eta[12]),
         x = expression(eta[11]),
         y = expression(eta[12])) +
    guides(color = guide_legend(title = expression("Function " ~ f(x[1]))),
           shape = guide_legend(title = expression("Function " ~ f(x[1])))) +
    scale_color_manual(labels = c(expression(lambda[1] ~ "          "),
                                  expression("  " ~ lambda[1] ~ frac(K, (K + x[1])))),
                       values = c("#fc9803", "#248aff")) +
    scale_shape_manual(labels = c(expression(lambda[1] ~ "          "),
                                  expression("  " ~ lambda[1] ~ frac(K, (K + x[1])))),
                       values = c(19, 17)) +
    geom_smooth(data = no_outlier_stats, method = "lm", se = FALSE)
  
  # Plot eta 12 and eta 23
  covar_plot <- ggplot(combined_stats, aes(x = eta_12, y = eta_23,
                                           color = func, shape = func)) +
    geom_point(size = 2.5, alpha = 0.5) +
    labs(title = expression("Normalized Covariances " ~ eta[12] ~ "and" ~ eta[23]),
         x = expression(eta[12]),
         y = expression(eta[23])) +
    guides(color = guide_legend(title = expression("Function " ~ f(x[1]))),
           shape = guide_legend(title = expression("Function " ~ f(x[1])))) +
    scale_color_manual(labels = c(expression(lambda[1] ~ "          "),
                                  expression("  " ~ lambda[1] ~ frac(K, (K + x[1])))),
                       values = c("#fc9803", "#248aff")) +
    scale_shape_manual(labels = c(expression(lambda[1] ~ "          "),
                                  expression("  " ~ lambda[1] ~ frac(K, (K + x[1])))),
                       values = c(19, 17)) +
    geom_smooth(data = no_outlier_stats, method = "lm", se = FALSE)
  
  # combine into subplots
  covar_subplots <- ggarrange(var_covar_plot, covar_plot, nrow = 1,
                              widths = c(1, 1), common.legend = TRUE,
                              legend = "right")
  
  return(covar_subplots)
}

# Remove outlier points
search_A_filtered <- search_A$numerical_stats[which(search_A$numerical_stats$eta_12 > 0),]
search_B_filtered <- search_B$numerical_stats[which(search_B$numerical_stats$eta_12 > 0),]

# Perform linear regression to find equation of line of best fit
lm(eta_23 ~ eta_12, data = search_A_filtered)$coefficients
lm(eta_23 ~ eta_12, data = search_B_filtered)$coefficients

# Create separate plots for analytical and numeric etas
analytical_subplots <- plotCovariance(points_A = search_A$analytical_stats,
                                      points_B = search_B$analytical_stats,
                                      no_outliers_A = search_A$analytical_stats,
                                      no_outliers_B = search_B$analytical_stats)
numeric_subplots <- plotCovariance(points_A = search_A$numerical_stats,
                                   points_B = search_B$numerical_stats,
                                   no_outliers_A = search_A_filtered,
                                   no_outliers_B = search_B_filtered)

# Save subplots
pdf("Analytical_Covariance.pdf", width = 13, height = 7)
annotate_figure(analytical_subplots, top = text_grob("Analytical Results", 
                                                     face = "bold", size = 14))
dev.off()

pdf("Numeric_Covariance.pdf", width = 13, height = 7)
annotate_figure(numeric_subplots, top = text_grob("Numerical Results", 
                                                  face = "bold", size = 14))
dev.off()