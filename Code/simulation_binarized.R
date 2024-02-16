# import packages and functions
source("packages.R")
source("ml_wrapper.R")
source("eval_functions.R")

set.seed(42)
options(scipen= 999)


# Set the number of simulations
n_rep <- 2
n_folds <- 4

# Set your list of classifiers
classifiers <- c(
  "adaboost",  "bagging", # "knn",
  "lda", "logit", "logit_nnet",
  "mlpc", "nb_bernulli", "nb_gaussian",
  "probability_forest", "qda", "ranger", "svm", "xgboost"
)

#datanames <- c("Imbens_2016","Linden_2015", "Archarky_2023")
datanames <- c("Archarky_2023")

# OvR Simulations
for (dataname in datanames){
  
  # Define the dimensions of the nested structure
  source(paste0("data/data_", dataname,".R"))
  
  n = nrow(X)
  t = length(unique(W))
  
  results <- vector("list", length = length(classifiers) + 1)
  names(results)[1:(length(classifiers) + 1)] <- c(classifiers, "outcome")
  
  # Loop through simulations and then through classifiers
  for (i in 1:n_rep) {
    print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
    
    # Define the dimensions of the nested structure
    source(paste0("data/data_", dataname,".R"))
    n = nrow(X)
    # store generated data set
    results[["outcome"]][[i]] <- dataset
    
    # Loop through classifiers
    for (c in seq_along(classifiers)) {
      classifier <- classifiers[c] 
      print(classifier)
      
      # define cross folding index
      fold = sample(1:n_folds,n,replace=T)
      predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
      colnames(predictions) <- seq(1,t,1)
      
      # form cross-fitted prediction for e
      for (f in 1:n_folds){
        print(paste(i,classifier,"fold ",f))
        # fit model 
        model <- do.call(ovr_fit, list(x = X[fold != f,], y = W[fold != f], 
                                       method = classifier 
                                       ))
        
        # predict
        predictions[fold == f,] <- do.call(predict.ovr_fit,
                                           list(model, x=X[fold != f,], y=W[fold != f], xnew = X[fold == f,],
                                                method = classifier))
        
      }
    ## Reorder the columns in the dataframe
    predictions <- predictions[, order(names(predictions))]
    
    # Store results
    results[[classifier]][[i]] <- predictions
    }
  }
  
  # store result as .rds file
  saveRDS(results, file = paste0("sim_results/simulation_results_ovr_", dataname,"_", n_rep, ".rds"))
  
}


start_time <- Sys.time()  # Record end time

# OvO simulations
for (dataname in datanames){
  
  # Define the dimensions of the nested structure
  source(paste0("data/data_", dataname,".R"))
  
  n = nrow(X)
  t = length(unique(W))
  
  results <- vector("list", length = length(classifiers) + 1)
  names(results)[1:(length(classifiers) + 1)] <- c(classifiers, "outcome")
  
  # Loop through simulations and then through classifiers
  for (i in 1:n_rep) {
    print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
    
    # Define the dimensions of the nested structure
    source(paste0("data/data_", dataname,".R"))
    n = nrow(X)
    # store generated data set
    results[["outcome"]][[i]] <- dataset
    
    # Loop through classifiers
    for (c in seq_along(classifiers)) {
      classifier <- classifiers[c] 
      print(classifier)
      
      # define cross folding index
      fold = sample(1:n_folds,n,replace=T)
      predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
      colnames(predictions) <- seq(1,t,1)
      
      # form cross-fitted prediction for e
      for (f in 1:n_folds){
        print(paste(i,classifier,"fold ",f))
        # fit model 
        model <- do.call(ovo_fit, list(x = X[fold != f,], y = W[fold != f], 
                                       method = classifier 
        ))
        
        # predict
        predictions[fold == f,] <- do.call(predict.ovo_fit,
                                           list(model, x=X[fold != f,], y=W[fold != f], xnew = X[fold == f,],
                                                method = classifier))
        
      }
      ## Reorder the columns in the dataframe
      predictions <- predictions[, order(names(predictions))]
      
      # Store results
      results[[classifier]][[i]] <- predictions
    }
  }
  
  # store result as .rds file
  saveRDS(results, file = paste0("sim_results/simulation_results_ovo_", dataname,"_", n_rep, ".rds"))
  
}


end_time <- Sys.time()  # Record end time
elapsed_time <- end_time - start_time
cat(n_rep, "OvO simulations:", "elapsed time:", elapsed_time, "\n")

