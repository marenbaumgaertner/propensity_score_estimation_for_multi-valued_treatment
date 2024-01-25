

brier_score <- function(probabilities, outcome, binary = TRUE){
  # reshape probabilities
  #if(binary==TRUE){
  #  if(!is.matrix(probabilities)){
  #    probabilities <- as.matrix(probabilities)
  #  }
  #  if(dim(probabilities)[2]>1){
  #    probabilities <- as_tibble(probabilities)
  #    if("1" %in% colnames(probabilities)){
  #      probabilities <- as.matrix(probabilities$`1`)
  #    }else if("V2" %in% colnames(probabilities)){
  #      probabilities <- as.matrix(probabilities$`V2`)
  #    }else{
  #      # @Maren: some prediction have V1 and V2 as colnames, add this to code.
  #      # @Maren: too much characer for 1 line, check out how you can do elegant line breaks
  #      paste("Probabilities argument does not indicates which column contains the probabilities for class 1. Rename column of just use corresponding column as input.")
  #      break
  #    }
  #  }
  #  
  #}#else{
  
  #}
  
  # reshape outcome for multi-valued treatment from 1xN to KxN
  if(length(unique(outcome))==2){
    outcome <- as.matrix(outcome)
  }else{
    outcome <- class.ind(as.matrix(outcome))
  }
  
  # @Maren: wie kann ich R sagen dass es okay ist wenn beide dimension identisch sind? -> unique :P
  if(unique(dim(probabilities) == dim(outcome))){
    bs <- mean(rowMeans((probabilities - outcome)^2, na.rm = TRUE), na.rm = TRUE)
  }
  return(bs)
  
}



make_calibtation_plot <- function(probabilities, outcome, method) {
  
  class_names = colnames(probabilities)
  
  if (length(unique(outcome))==2){
    calibration_data <- cbind(probabilities, outcome)
    colnames(calibration_data) <- c("predicted_probabilities", "actual_outcomes")
    
  }else{
    probabilities <- probabilities %>% 
      pivot_longer(cols = class_names, 
                   names_to = "group",
                   values_to = "predicted_probabilities")
    true_outcome <- class.ind(W_test) %>% 
      as_tibble() %>% 
      pivot_longer(cols = class_names, 
                   names_to = "group",
                   values_to = "actual_outcomes")
    calibration_data <- cbind(probabilities, actual_outcomes=true_outcome$actual_outcomes)
    
  }
  
  # Create probability bins (e.g., 10 bins)
  n_bins <- 10
  bin_width <- 1 / n_bins
  bin_centers <- seq(0, 1, by = bin_width) + bin_width / 2
  
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
  if (length(unique(W_test))==2){
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
      theme_bw()
    
  }else{
    custom_palette <- c("#8E063B", "#79D359", "#4B0055", "#F7A72B", "#FFBEC1")
    plot <- ggplot(calibration_summary, aes(x = mean_predicted_prob, y = mean_actual_outcome, 
                                            color=group)) +
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
      scale_color_manual(values = custom_palette) +
      theme_bw()
    
  }
  
  return(plot)
}


finetuning_soft_classifier <- function(x, y, method, iter = 10, n_folds=5){
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
    params <- list()
    # Store tuning parameters
    for (col_name in colnames(method$grid)) {
      if (is.factor(method$grid[[col_name]])) {
        params[[col_name]] <- as.character(method$grid[[col_name]][i])
      } else {
        params[[col_name]] <- method$grid[[col_name]][i]
      }
    }
    bs_iter <- array(0, dim=iter)
    for (j in 1:iter) {
      
      n = nrow(dataset)
      
      # define cross folding index
      fold = sample(1:n_folds,n,replace=T)
      predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
      
      # cross validation
      for (f in 1:n_folds){
        
        # fit model 
        model <- do.call(method$fit, list(x = X[fold != f,], y = W[fold != f]))
        
        # predict
        predictions[fold == f,] <- do.call(method$predict,
                                           list(model, X[fold != f,], W[fold != f], xnew = X[fold == f,]))
        
      }
      bs_iter[j] <- brier_score(probabilities = predictions, outcome = y)
      
      
    }
    # only store average over iterations
    bs_mat[i, 'BS'] <- mean(bs_iter)
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
