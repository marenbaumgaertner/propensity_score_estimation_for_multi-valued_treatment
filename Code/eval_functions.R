

brier_score <- function(probabilities, outcome, binary = TRUE){

  # reshape outcome for multi-valued treatment from 1xN to KxN
  if(length(unique(outcome))==2){
    outcome <- as.matrix(outcome)
  }else{
    outcome <- class.ind(as.matrix(outcome))
  }
  
  if(unique(dim(probabilities) == dim(outcome))){
    bs <- mean(rowMeans((probabilities - outcome)^2, na.rm = TRUE), na.rm = TRUE)
  }

  return(bs)
  
}

cross_entropy <- function(probabilities, outcome) {
  #  Adjust indices since R indexing starts from 1 if necessary
  if (min(outcome)==0) outcome = outcome + 1
  
  if(unique(ncol(probabilities) == length(unique(outcome)))){
    # Select the predicted probabilities corresponding to true classes
    probs <- probabilities[cbind(1:n, outcome)] 
    
    # compute logarithm loss and divide by number of observations
    ll = -sum(log(probs+1e-20))/length(outcome)
  }else{
    stop("Error: Probabilities matrix does not match outcome classes.")
  }
  
  
  return(ll)
}


make_calibtation_plot <- function(probabilities, outcome, method, n_bins = 10) {
  
  class_names = colnames(probabilities)
  
  if (length(unique(outcome))>2) if (min(outcome)==0) outcome = outcome+1
  
  if (length(unique(outcome))==2){
    calibration_data <- cbind(probabilities, outcome)
    colnames(calibration_data) <- c("predicted_probabilities", "actual_outcomes")
    
  }else{
    probabilities <- probabilities %>% 
      pivot_longer(cols = class_names, 
                   names_to = "group",
                   values_to = "predicted_probabilities")
    true_outcome <- class.ind(outcome) %>% 
      as_tibble() %>% 
      pivot_longer(cols = class_names, 
                   names_to = "group",
                   values_to = "actual_outcomes")
    calibration_data <- cbind(probabilities, actual_outcomes=true_outcome$actual_outcomes)
    
  }
  
  # Create probability bins (e.g., 10 bins)
  bin_width <- 1 / n_bins
  bin_centers <- seq(0, 1, by = bin_width) + bin_width / 2
  bin_centers <- bin_centers[1:length(bin_centers)-1]
  
  # Group data into bins
  calibration_data$probability_bin <- cut(
    calibration_data$predicted_probabilities,
    breaks = c(-Inf, bin_centers, Inf), # use Inf such that very small/large probabilities get considered
    labels = FALSE
  )  
  
  # Calculate mean predicted probability and mean actual outcome for each bin
  calibration_summary <- calibration_data %>%
    dplyr::group_by(probability_bin, group) %>%
    dplyr::summarize(
      mean_predicted_prob = mean(predicted_probabilities),
      mean_actual_outcome = mean(actual_outcomes)
    )
  
  # Create a calibration plot
  if (length(unique(outcome))==2){
    plot <- ggplot(calibration_summary, aes(x = mean_predicted_prob, y = mean_actual_outcome)) +
      geom_point() +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Identity line
      labs(
        x = "Mean Predicted Probability",
        y = "Mean Actual Outcome",
        title = paste0("Calibration Plot of ", method)
      ) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_color_paletteer_d("palettetown::delcatty", direction = 1, dynamic = FALSE) +
      theme_bw()
    
  }else{
    custom_palette <- c("#D23736", "#79D359", "#4B0055", "#F7A72B", "#FFBEC1", 
                        "#8E063B", "#374f3f", "#631b1f", "#c45400", "#842c5c",
                        "#910e08", "#1f3611", "#6c2132")
    plot <- ggplot(calibration_summary, aes(x = mean_predicted_prob, y = mean_actual_outcome, 
                                            color=group)) +
      geom_point() +
      geom_line() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Identity line
      labs(
        x = "Mean Predicted Probability",
        y = "Mean Actual Outcome",
        color = "Treatment",
        title = paste0("Calibration Plot of ", method)
      ) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      scale_color_manual(values = custom_palette) +
      theme_bw()
    
  }
  
  return(plot)
}

kl_convergence <- function(p, q_matrix) {
  epsilon <- 1e-10
  p <- p / sum(p) # Normalize p to sum to 1
  K <- length(p)
  total_divergence <- 0

  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      u_ij <- p[i] / (p[i] + p[j] + epsilon)
      r_ij <- q_matrix[i, j]
      
      if (!is.na(r_ij)) {
        total_divergence <- total_divergence + 
          r_ij * log((r_ij + epsilon) / (u_ij + epsilon)) + 
          (1 - r_ij) * log(((1 - r_ij) + epsilon) / ((1 - u_ij) + epsilon)) 
      }
    }
  }
  
  return(total_divergence)
}


finetuning_soft_classifier <- function(x, y, method, n_folds=5){
  #library(method$library) # install and load package
  
  if (is.null(method$grid)) {
    stop("Error: No grid defined")
  }
  if (is.null(method$fit)) {
    stop("Error: No fit function defined")
  }
  if (is.null(method$predict)) {
    stop("Error: No predict function defined")
  }
  
  set.seed(42)
  grid <- method$grid
  grids <- nrow(grid)
  
  # initialize empty matrix for resulting Brier Score
  bs_mat <- grid
  bs_mat$BS <- NA
  
  # bind input data
  dataset <- data.frame(y = y, x = x)
  
  for (i in 1:grids) {
    print(paste0("grid ",i, " of ", grids))
    params <- list()
    # Store tuning parameters
    for (col_name in colnames(method$grid)) {
      if (is.factor(method$grid[[col_name]])) {
        params[[col_name]] <- as.character(method$grid[[col_name]][i])
      } else {
        params[[col_name]] <- method$grid[[col_name]][i]
      }
    }
      
    n = nrow(dataset)
    
    # define cross folding index
    fold = sample(1:n_folds,n,replace=T)
    predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
    
    # cross validation
    for (f in 1:n_folds){
      # fit model 
      model <- do.call(method$fit, list(x = X[fold != f,], y = W[fold != f], params))
      #print(paste0("Model param lambda = ", model$lambda))
      #print(paste0("Grid param lambda = ", params$lambda))
      
      # predict
      predictions[fold == f,] <- do.call(method$predict,
                                         list(model, X[fold != f,], W[fold != f], xnew = X[fold == f,]))
      
    }
      
    # only store brier score
    bs_mat[i, 'BS'] <- brier_score(probabilities = predictions, outcome = y)
  }
  
  best_params <- bs_mat[which.min(bs_mat$BS), ]
  
  
  list('performance_mat' = bs_mat,
       'best_params' = best_params)
}

#true_outcome <- class.ind(W_test) %>% as_tibble()
#probabilities <- results$bagging$predictions
#class_names = colnames(probabilities)


#probabilities <- probabilities %>% 
#  pivot_longer(cols = class_names, 
#               names_to = "group",
#               values_to = "predicted_probabilities")
#true_outcome <- class.ind(W_test) %>% 
#  as_tibble() %>% 
#  pivot_longer(cols = class_names, 
#               names_to = "group",
#               values_to = "actual_outcomes")

#calibration_data <- cbind(probabilities, actual_outcomes=true_outcome$actual_outcomes)

## Create probability bins (e.g., 10 bins)
#n_bins <- 10
#bin_width <- 1 / n_bins
#bin_centers <- seq(0, 1, by = bin_width) + bin_width / 2

## Group data into bins
## Group data into bins
#calibration_data$probability_bin <- cut(
#  calibration_data$predicted_probabilities,
#  breaks = c(-Inf, bin_centers, Inf),
#  labels = FALSE
#)
## Calculate mean predicted probability and mean actual outcome for each bin
#calibration_summary <- calibration_data %>%
#  dplyr::group_by(probability_bin, group) %>%
#  dplyr::summarize(
#    mean_predicted_prob = mean(predicted_probabilities),
#    mean_actual_outcome = mean(actual_outcomes)
#  )


## Create a calibration plot
#plot <- ggplot(calibration_summary, aes(x = mean_predicted_prob, y = mean_actual_outcome, color=group)) +
#  geom_point() +
#  geom_line() +
#  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Identity line
#  labs(
#    x = "Mean Predicted Probability",
#    y = "Mean Actual Outcome",
#    title = paste0("Calibration Plot of ", method)
#  ) +
#  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#  theme_bw()
#plot#
