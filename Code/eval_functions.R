

brier_score <- function(probabilities, outcome, binary = TRUE){
  # reshape probabilities
  if(binary==TRUE){
    if(!is.matrix(probabilities)){
      probabilities <- as.matrix(probabilities)
    }
    if(dim(probabilities)[2]>1){
      probabilities <- as_tibble(probabilities)
      if("1" %in% colnames(probabilities)){
        probabilities <- as.matrix(probabilities$`1`)
      }else if("V2" %in% colnames(probabilities)){
        probabilities <- as.matrix(probabilities$`V2`)
      }else{
        # @Maren: some prediction have V1 and V2 as colnames, add this to code.
        # @Maren: too much characer for 1 line, check out how you can do elegant line breaks
        paste("Probabilities argument does not indicates which column contains the probabilities for class 1. Rename column of just use corresponding column as input.")
        break
      }
    }
    
  }#else{
  
  #}
  
  # reshape outcome for multi-valued treatment from 1xN to KxN
  if(binary == TRUE){
    outcome <- as.matrix(outcome)
  }else{
    outcome <- class.ind(as.matrix(outcome))
  }
  
  # @Maren: wie kann ich R sagen dass es okay ist wenn beide dimension identisch sind?
  if(unique(dim(probabilities) == dim(outcome))){
    bs <- mean((probabilities - outcome)^2)
    
  }
  return(bs)
  
}



make_calibtation_plot <- function(probabilities, true_outcome, method) {
  ## Calibration plot
  #if (ncol(probabilities)>)
  probabilities <- as_tibble(probabilities)
  calibration_data <- cbind(probabilities, true_outcome)
  colnames(calibration_data) <- c("predicted_probabilities", "actual_outcomes")
  
  # Create probability bins (e.g., 10 bins)
  n_bins <- 10
  bin_width <- 1 / n_bins
  bin_centers <- seq(0, 1, by = bin_width) + bin_width / 2
  
  # Group data into bins
  calibration_data$probability_bin <- cut(calibration_data$predicted_probabilities, breaks = bin_centers, labels = FALSE)
  
  # Calculate mean predicted probability and mean actual outcome for each bin
  calibration_summary <- calibration_data %>%
    dplyr::group_by(probability_bin) %>%
    dplyr::summarize(
      mean_predicted_prob = mean(predicted_probabilities),
      mean_actual_outcome = mean(actual_outcomes)
    )
  
  # Create a calibration plot
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
  
  return(plot)
}
