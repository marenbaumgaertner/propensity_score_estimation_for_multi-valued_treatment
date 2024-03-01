# import packages and functions
source("packages.R")
source("ml_wrapper.R")
source("eval_functions.R")

set.seed(42)
options(scipen= 999)


# Set the number of simulations
n_rep <- 100
n_folds <- 4

# Set your list of classifiers
classifiers <- c(
  "adaboost",  "bagging", "knn",
  "lda", "logit", "logit_nnet",
  "mlpc", "nb_bernulli", "nb_gaussian",
  "probability_forest", "qda", "ranger", "svm", "xgboost"
)

# define datasets
datanames <- c("Imbens_2016", "Linden_2015", "Acharky_2023")



# Add OvR Simulations to  multi-class simulation results
for (dataname in datanames){
  
  # Define the dimensions of the nested structure
  results <- readRDS(paste0("sim_results/simulation_results_", dataname,"_", n_rep, ".rds"))
 
  
  # Loop through simulations and then through classifiers
  for (i in 1:n_rep) {
    print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
    
    # Extract data from list
    dataset <- results[["outcome"]][[i]]
    W <- dataset$W
    X <- dataset %>% select(!c(W,Y)) %>%  as.matrix()
    
    n = nrow(X)
    t = length(unique(W))

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
        
        tryCatch({
          predictions[fold == f,] <- do.call(predict.ovr_fit,
                                             list(model, x=X[fold != f,], y=W[fold != f], xnew = X[fold == f,],
                                                  method = classifier))
        }, error = function(e) {
          # Handle the error here
          print(paste("An error occurred:", e))
          predictions[fold == f,] <- rep(NA, length(unique(W)))
        })
        
      }
      ## Reorder the columns in the dataframe
      predictions <- predictions[, order(names(predictions))]
      
      # Store results
      results[[paste0("ovr_", classifier)]][[i]] <- predictions
    }
  }
  # store result as .rds file
  saveRDS(results, file = paste0("sim_results/simulation_results_ovr_", dataname,"_", n_rep, ".rds"))
}

# Add OvO Simulations to OvR and  multi-class simulation results
datanames <- c("Imbens_2016", "Linden_2015")

for (dataname in datanames){
  
  # Load result list
  results <- readRDS(paste0("sim_results/simulation_results_ovr_", dataname,"_", 100, ".rds"))
  
  
  # Loop through simulations and then through classifiers
  for (i in 1:n_rep) {
    print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
    
    # Extract data from list
    dataset <- results[["outcome"]][[i]]
    W <- dataset$W
    X <- dataset %>% select(!c(W,Y)) %>%  as.matrix()
    
    n = nrow(X)
    t = length(unique(W))

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
        
        tryCatch({
          predictions[fold == f,] <- do.call(predict.ovo_fit,
                                             list(model, x=X[fold != f,], y=W[fold != f], xnew = X[fold == f,],
                                                  method = classifier))
        }, error = function(e) {
          # Handle the error here
          print(paste("An error occurred:", e))
          predictions[fold == f,] <- rep(NA, length(unique(W)))
        })
        
      }
      ## Reorder the columns in the dataframe
      predictions <- predictions[, order(names(predictions))]
      
      # Store results
      results[[paste0("ovo_", classifier)]][[i]] <- predictions
    }
  }
  # store result as .rds file
  saveRDS(results, file = paste0("sim_results/simulation_results_ovo_", dataname,"_", n_rep, ".rds"))
}



# Execution time sequential  1359.77       10.62     1570.08    26 Minuten
# 2546.28       17.67     2604.61 43 Minuten

# Parallelize computations for Acharky
datanames <- c("Acharky_2023")


system.time({
  for (dataname in datanames){
    
    # Define the dimensions of the nested structure
    results <- readRDS(paste0("sim_results/simulation_results_ovo_", dataname,"_", 100, ".rds"))
    
    # Loop through simulations and then through classifiers
    for (i in 1:n_rep) {
      print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
      
      # Extract data from list
      dataset <- results[["outcome"]][[i]]
      W <- dataset$W
      X <- dataset %>% select(!c(W,Y)) %>%  as.matrix()
      
      n = nrow(X)
      t = length(unique(W))
      
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
          
          cl <- makeCluster(8)
          registerDoParallel(cl)
          clusterExport(cl, c(paste0(classifier, "_fit"), 
                              paste0("predict.",classifier, "_fit"),
                              "glmnet", "nnet", "ranger", 
                              "probability_forest", "multinom", "train.kknn",
                              "gaussian_naive_bayes", "bernoulli_naive_bayes",
                              "adaboost", "bagging", "svm", "lda", "qda", "xgboost",
                              "class.ind", "naive_bayes", "rpart.control",
                              "%>%", "kl_convergence", "optim"
                              
          ))
          
          # fit model 
          model <- do.call(ovo_fit_parallel, list(x = X[fold != f,], y = W[fold != f], 
                                         method = classifier 
          ))
          
          
          #stop()
          tryCatch({
            predictions[fold == f,] <- do.call(predict.ovo_fit_parallel,
                                               list(model, x=X[fold != f,], y=W[fold != f], xnew = X[fold == f,],
                                                    method = classifier))
          }, error = function(e) {
            # Handle the error here
            print(paste("An error occurred:", e))
            predictions[fold == f,] <- rep(NA, length(unique(W)))
          })
          stopCluster(cl)
        }
        ## Reorder the columns in the dataframe
        predictions <- predictions[, order(names(predictions))]
        
        # Store results
        results[[paste0("ovo_", classifier)]][[i]] <- predictions
      }
    }
    # store result as .rds file
    saveRDS(results, file = paste0("sim_results/simulation_results_ovo_", dataname,"_", n_rep, ".rds"))
  }
})
