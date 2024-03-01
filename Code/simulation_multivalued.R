# import packages and functions
source("packages.R")
source("ml_wrapper.R")
source("eval_functions.R")

set.seed(42)
options(scipen= 999)


# Set the number of simulations and number of cross folds
n_rep <- 100
n_folds <- 5

# Set your list of classifiers
classifiers <- c(#"ovo", "ovr",
  "knn", "lda", "logit", "logit_nnet",
  "mlpc", "nb_bernulli", "nb_gaussian", 
  "probability_forest", "ranger", "qda" #"xgboost"
)

#classifiers <- c("ranger")

#dataname <- "Imbens_2016" ## 1 Sim 15 s
#dataname <- "Linden_2015" ## 1 sim 10 s
#dataname <- "Archarky_2023"

datanames <- c("Imbens_2016","Linden_2015", "Acharky_2023")
#datanames <- c("Imbens_2016","Linden_2015")

for (dataname in datanames){
  
 
  tuning_results <- readRDS(paste0("sim_results/final_tuning_results_", dataname, ".rds"))
  tuning_results$xgboost = tuning_results$xgboost_method
  tuned_classifiers <- names(tuning_results)
  
  
  # Define the dimensions of the nested structure
  source(paste0("data/data_", dataname,".R"))
  
  n = nrow(X)
  t = length(unique(W))
  
  results <- vector("list", length = length(classifiers) + 1)
  names(results)[1:(length(classifiers) + 1)] <- c(classifiers, "outcome")
  
  # track overall execution time
  start_time <- Sys.time()  # Record start time
  
  # Loop through simulations and then through classifiers
  for (i in 1:n_rep) {
    print(paste0("iteration ",i, " of ", n_rep)) # to keep track of execution
    
    # generate new data for each iteration
    source(paste0("data/data_", dataname,".R"))
    n = nrow(X)
    print(head(dataset))
    
    # store generated data set
    results[["outcome"]][[i]] <- dataset
  
    for (c in seq_along(classifiers)) {
      classifier <- classifiers[c]
      #print(classifier)
    
      if (classifier %in% tuned_classifiers){
        params_raw <- tuning_results[[classifier]]
        params_raw$BS <- NULL
        rownames(params_raw) <- NULL
        params <- list()
        # Store tuning parameters
        for (col_name in colnames(params_raw)) {
          if (is.factor(params_raw[[col_name]])) {
            params[[col_name]] <- as.character(params_raw[[col_name]][1])
          } else {
            params[[col_name]] <- params_raw[[col_name]][1]
          }
        }
      }else if(classifier %in% c("ovo","ovr")){
        params = list(method="logit")
      }else{
        params = list()
      }
      
      # define cross folding index
      fold = sample(1:n_folds,n,replace=T)
      predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
      colnames(predictions) <- seq(1,t,1)
      
      # form cross-fitted prediction for e
      for (f in 1:n_folds){
        
        # fit model 
        model <- do.call(paste0(classifier, "_fit"), list(x = X[fold != f,], y = W[fold != f], params))
        
        # predict
        predictions[fold == f,] <- do.call(paste0("predict.", classifier, "_fit"),
                                           list(model, X[fold != f,], W[fold != f], xnew = X[fold == f,]))
        
      }
      # Store the results
      results[[classifier]][[i]] <- predictions
    }
   }
  
  end_time <- Sys.time()  # Record end time
  elapsed_time <- end_time - start_time
  cat(n_rep, "Simulations:", "elapsed time:", elapsed_time, "\n")
  
  # save the results as an RDS file
  saveRDS(results, file = paste0("sim_results/simulation_results_", dataname,"_", n_rep, ".rds"))

}





## Loop through classifiers and then through iteration -> dumm
#for (c in seq_along(classifiers)) {
#  classifier <- classifiers[c]
#  print(classifier)
#  
#  if (classifier %in% tuned_classifiers){
#    params_raw <- tuning_results[[classifier]]$best_params
#    params_raw$BS <- NULL
#    rownames(params_raw) <- NULL
#    params <- list()
#    # Store tuning parameters
#    for (col_name in colnames(params_raw)) {
#      if (is.factor(params_raw[[col_name]])) {
#        params[[col_name]] <- as.character(params_raw[[col_name]][1])
#      } else {
#        params[[col_name]] <- params_raw[[col_name]][1]
#      }
#    }
#  }else{
#    params = list()
#  }
#  
#  # Initialize a list to store results for each simulation
#  results_list <- vector("list", length = n_rep)
#  
#  start_time_c <- Sys.time()  # Record start time
#  
#  # Loop through simulations
#  for (i in 1:n_rep) {
#    # use data from Imbens (2015)
#    source(paste0("data/data_", dataname,".R"))
#    
#    n = nrow(X)
#    
#    # define cross folding index
#    fold = sample(1:n_folds,n,replace=T)
#    predictions = data.frame(matrix(NA,nrow=n, ncol = length(unique(W))))
#    
#    # form cross-fitted prediction for e
#    for (f in 1:n_folds){
#
#      # fit model 
#      model <- do.call(paste0(classifier, "_fit"), list(x = X[fold != f,], y = W[fold != f], params))
#      
#      # predict
#      predictions[fold == f,] <- do.call(paste0("predict.", classifier, "_fit"),
#                                        list(model, X[fold != f,], W[fold != f], xnew = X[fold == f,]))
#      
#    }
#  
#    # Store the results
#    results_list[[i]] <- predictions
#  }
#  end_time_c <- Sys.time()  # Record end time
#  elapsed_time_c <- end_time_c - start_time_c
#  print(paste("Classifier", classifiers[c], "elapsed time:", elapsed_time_c))
#  
#  # Store the results for the current classifier
#  names(all_results)[c] <- classifier
#  all_results[[c]] <- results_list
#}
