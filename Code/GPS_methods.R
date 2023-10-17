# Import libraries and dataset
library(dplyr)
library(reshape2)
library(randomForest)
library(ggplot2)
library(caTools)
library(mlogit)
library(nnet)
library(caret)
library(glmnet)
library(naivebayes)
library(rpart)
library(randomForest)
library(MASS)
library(FNN)
library(e1071)
library(igraph)
library(caret)
library(sparklyr)
library(RSSL)
library(sgd)
library(causalDML)
#source("MetaLearners_tools.R");rm(cl, nc)
source("eval_functions.R")
set.seed(42)


# create toy dataset (linear DGP)
outcome <- function(t, X){
  return( (1+t)*X )
}

n <- 5000
K <- 5

X1 <- runif(n = n)
X2 <- runif(n = n)
X3 <- runif(n = n)
W <- sample(seq(1, K, by = 1), replace = TRUE, n)

W.levels <- sort(unique(W))
K == length(W.levels)

#Y <- outcome(W, X) +  rnorm(n, sd = 0.05)
#dataset <- data.frame(X, W = W, Y = Y)
dataset <- data.frame(X1, X2, X3, W = W)

split_indices <- sample.split(dataset$W, SplitRatio = 0.8)
training_data <- dataset[split_indices, ]
test_data <- dataset[!split_indices, ]
W_training <- data.frame(W = training_data$W)
X_training <- data.frame(training_data[ ,1:3])
W_test <- data.frame(W = test_data$W)
X_test <- data.frame(test_data[,1:3])

formula <- W ~ X1 + X2


K_GPS_unif <- function(x){
  W = x[1]
  X = x[2]
  k = W*(K-1)
  if( X > qunif(k/K) & X <= qunif((k+1)/K)){
    GPS = (K+1)/(2*K)
  }else{
    GPS = 1/(2*K)
  }
  return(GPS)
}


# 0. uniformly distributed propensit score
# @Maren: das tut nicht, to be corrected later
e_hat_unif <- matrix(NA, nrow = nrow(test_data), ncol = K)
for(k in 1:K){
  e_hat_unif[,k] = apply(data.frame(W.levels[k], test_data$X1), 1, K_GPS_unif)
}

# 1. logistic regression
#unpenalized mlogit via glmnet is far faster than the mlogit package
e_hat_logit <- matrix(NA, nrow = nrow(test_data), ncol = K)
# model_logit = mlogit::mlogit(W ~ X, data = training_data) # need the data to be reshaped 
model_logit <- glmnet::mlogit(formula, data = training_data, model = TRUE) # trace=FALSE?
model_logit <- glm(formula, data = training_data, model = TRUE, family = "multinomial")
model_logit <- nnet::multinom(formula, data = training_data, model = TRUE)
e_hat_logit <- predict(model_logit, newdata = X_test, type = "prob")

bs_logit <- brier_score(e_hat_logit, W_test, binary = FALSE)

# 2. logistic regression with cross-validation
e_hat_logit_cv <- matrix(NA, nrow = nrow(test_data), ncol = K)
ctrl <- trainControl(
  method = "cv",  # Cross-validation method
  number = 5      # Number of folds
)

#model_logit_cv <- train(
#  formula,
#  data = training_data,
#  method = "multinom",  # Specify the method
#  trControl = ctrl     # Use the defined cross-validation settings
#)
#
#e_hat_logit_cv <- predict(model_logit_cv, newdata = X_test, type = "prob")

# 3. Ridge classification
model_ridge <- glmnet(x = as.matrix(X_training), y = as.matrix(W_training), alpha = 0, family = "multinomial")
e_hat_ridge <- predict(model_ridge, newx = as.matrix(X_test), s = 0.01, type = "response")

# 4. Ridge classification cross-validation
model_ridge_cv <- cv.glmnet(x = as.matrix(X_training), y = as.matrix(W_training), alpha = 0, family = "multinomial")
e_hat_ridge_cv <- predict(model_ridge_cv, newx = as.matrix(X_test), s = 0.01, type = "response")

# 5. Naive bayes classifier (Bernulli)
model_nb_bernulli <- bernoulli_naive_bayes(X_training, as.factor(training_data$W))
e_hat_nb_bernulli <- predict(model_nb_bernulli, newdata = as.matrix(X_test), type = "prob")

# 6. Naive bayes classifier (Gaussian)
model_nb_gaussian <- naive_bayes(X_training, as.factor(training_data$W))
e_hat_nb_gaussian <- predict(model_nb_gaussian, newdata = as.matrix(X_test), type = "prob")

all(e_hat_nb_bernulli == e_hat_nb_gaussian)

# 7. Decision Tree classifier
model_decision_tree <- rpart(formula, data = training_data, method = "class")
e_hat_decision_tree <- predict(model_decision_tree, X_test, type = "prob")

# 8. Extra Tree classifier
# R package not available, just for binary classification
#model_extra_tree <-
#e_hat_extra_tree <-

# 8. XGBoost
model_xgb <- xgboost(data = as.matrix(training_data), label = as.matrix(training_data$W),
                     nrounds = 30, max_depth = 5, num_class = K,
                     eta = 0.1,
                     objective = "multi:softprob", eval_metric = "mlogloss")
e_hat_xgb_raw <- predict(model_xgb, as.matrix(test_data))

# check whether first 5 element sum uü to 1
e_hat_xgb_raw[1] + e_hat_xgb_raw[2] + e_hat_xgb_raw[3] + e_hat_xgb_raw[4] + e_hat_xgb_raw[5]

e_hat_xgb <- matrix(e_hat_xgb_raw, nrow = 1000, ncol = 5, byrow = TRUE)

brier_score_50 <- mean((e_hat_xgb-class.ind(test_data$W))^2)

# compute Brier Score for multivalued treatment
brier_score_10 <- mean((e_hat_xgb-class.ind(test_data$W))^2)

probabilities <- e_hat_rf
outcome<-W_test



# 9. Random Forest
e_hat_rf <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_rf <- randomForest(y = as.factor(training_data$W), x = X_training)
e_hat_rf <- predict(model_rf, X_test, type = "prob")

# 10. Label Propagation
#model_label_prop <- clustee_label_prop()
#e_hat_label_prop <-
  
GRFClassifier(as.matrix(X_training), as.factor(W_training), as.matrix(X_test))

## Perform Principal Component Analysis (PCA) to reduce dimensionality
#pca_result <- prcomp(X_training, scale. = TRUE)
#
## Extract the first few principal components as features
#num_components <- 2  # You can adjust the number of components
#reduced_features <- pca_result$x[, 1:num_components]
#
## Calculate the similarity matrix based on the reduced features
#similarity_matrix <- as.matrix(dist(reduced_features))# Create a graph from the similarity matrix
#graph <- graph.adjacency(similarity_matrix, mode = "undirected")
#
## Perform label propagation
#labels <- as.factor(W_training)  # W is the target variable
#
## Set the number of iterations (you can adjust this parameter)
#num_iterations <- 10
#
## Use label propagation algorithm
#for (i in 1:num_iterations) {
#  labels <- label_propagation(graph, labels)
#}
#
## Calculate probability estimates for each class
#probability_estimates <- table(labels) / length(labels)

# 11. Label Spreading
model_label_spread <-
e_hat_label_spread <-

# 12. Linear Support Vector Machine
e_hat_svm <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_svm <- svm(y = as.factor(training_data$W), x = X_training, probability = TRUE)
e_hat_svm <- predict(model_svm, X_test, probability = TRUE) %>% attr("probabilities")
bs_svm <- brier_score(e_hat_svm, W_test, binary = FALSE)

# 13. K-Nearest-Neighbors
e_hat_knn <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_knn <- FNN::knn(train = as.matrix(X_training), test = as.matrix(X_test), cl = training_data$W, prob = TRUE)
model_knn <- class::knn(train = as.matrix(X_training), test = as.matrix(X_test), cl = training_data$W, prob = TRUE)
e_hat_knn <- attributes(model_knn)$prob

# 14. Nearest Centroid
#model_nearest_centriod <-
#e_hat_nearest_centriod <-

# 15. Nearest-Neighbor (Radius)
#model_nn_radius <-
#e_hat_nn_radius <-

# 16. Linear Discriminant Analysis
model_lda <- lda(formula, data = training_data)
e_hat_lda <- predict(model_lda, X_test, type = "prob")$posterior
bs_lda <- brier_score(e_hat_lda, W_test, binary = FALSE)

# 17. Quadratic Discriminant Analysis
model_qda <- qda(formula, data = training_data)
e_hat_qda <- predict(model_qda, X_test, type = "prob")$posterior
bs_qda <- brier_score(e_hat_qda, W_test, binary = FALSE)


# 18. Multi-layer Perceptron classifier
model_mlpc <- nnet(x = scale(X_training), y = class.ind(training_data$W), linout = FALSE, size = c(10, 5), softmax = TRUE, maxit = 200)
e_hat_mlpc <- predict(model_mlpc, X_test)

# 19. Stochastic Gradient Descent (SGD)
model_sgd <- sgd(as.matrix(X_training), W_training, 
                 model = "lm",
                 #loss = "logistic",  # Use logistic loss for classification
                 #eta0 = 0.01,         # Learning rate (you can adjust this)
                 #nepochs = 100        # Number of epochs (you can adjust this)
                 )

e_hat_sgd <- predict_all(model_sgd, as.matrix(X_test), type = "prop") # too much cols
e_hat_sgd <- predict(model_sgd, as.matrix(X_test, type = "response")) # just 1 col


# Function to perform stochastic gradient descent for multiclass classification
# Inputs:
#   - X_train: Training feature matrix (n_samples x n_features)
#   - y_train: Training labels (n_samples)
#   - X_test: Test feature matrix (n_samples x n_features)
#   - learning_rate: Learning rate for SGD
#   - num_epochs: Number of epochs (iterations) for training
#   - num_classes: Number of classes in the classification problem
# Outputs:
#   - prob_estimates: Probability estimates for each class on the test data (n_samples x num_classes)
#   - weights: Learned weights for the model

multiclass_sgd <- function(X_train, y_train, X_test, learning_rate, num_epochs, num_classes) {
  
  # Initialize weights with zeros
  weights <- matrix(0, ncol = num_classes, nrow = ncol(X_train))
  
  # Number of training samples
  n_samples <- nrow(X_train)
  
  # Convert labels to one-hot encoding
  y_train_onehot <- matrix(0, nrow = n_samples, ncol = num_classes)
  for (i in 1:n_samples) {
    y_train_onehot[i, y_train[i]] <- 1
  }
  
  # Perform SGD for the specified number of epochs
  for (epoch in 1:num_epochs) {
    
    # Shuffle the training data for each epoch
    index <- sample(1:n_samples)
    X_train <- X_train[index, , drop = FALSE]
    y_train_onehot <- y_train_onehot[index, , drop = FALSE]
    
    # SGD for each training sample
    for (i in 1:n_samples) {
      xi <- X_train[i, , drop = FALSE]
      yi <- y_train_onehot[i, , drop = FALSE]
      
      # Compute softmax scores
      scores <- exp(xi %*% weights)
      probs <- scores / sum(scores)
      
      # Compute the gradient
      gradient <- t(xi) %*% (probs - yi)
      
      # Update weights using gradient descent
      weights <- weights - learning_rate * gradient
    }
  }
  
  # Calculate probability estimates for the test data
  scores_test <- exp(X_test %*% weights)
  prob_estimates <- scores_test / rowSums(scores_test)
  
  return(list(prob_estimates = prob_estimates, weights = weights))
}

# Example usage:
# Replace X_train, y_train, X_test, learning_rate, num_epochs, and num_classes with your data and hyperparameters.
# result <- multiclass_sgd(X_train, y_train, X_test, learning_rate = 0.01, num_epochs = 100, num_classes = 3)
# prob_estimates <- result$prob_estimates
# weights <- result$weights

e_hat_sgd <- multiclass_sgd(as.matrix(X_training), as.factor(W_training), as.matrix(X_test),
                            learning_rate = 0.01, num_epochs = 100, num_classes = K)$prob_estimates





library(caret)

# Function to create a one-vs-one classifier
fit_ovo_classifier <- function(data, labels, classifier = "svm") {
  class_labels <- unique(labels)
  num_classes <- length(class_labels)
  ovo_classifiers <- list()
  
  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes - 1)) {
    for (j in (i + 1):num_classes) {
      class_i <- class_labels[i]
      class_j <- class_labels[j]
      subset_data <- data[labels %in% c(class_i, class_j), ]
      subset_labels <- labels[labels %in% c(class_i, class_j)]
      
      # Train a binary classifier on the subset of data
      if (classifier == "svm") {
       # ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- train(
       #   x = subset_data,
       #   y = factor(subset_labels, levels = c(class_i, class_j)),
       #   method = "svmLinear",
       #   trControl = trainControl(method = "cv")
       # )
        # @Maren: change to ml_methods[[ml]]
        ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- svm(
          y = factor(subset_labels, levels = c(class_i, class_j)),
          x = subset_data,
          probability = TRUE)
      }
      # You can add other binary classifiers as needed
      
    }
  }
  
  return(ovo_classifiers)
}

ovo_classifier_list <- ovo_classifiers
newx <- X_test

predict_ovo_classifier <- function(ovo_classifier_list, newx, classifier = "svm") {
  n_classifiers <- length(ovo_classifier_list)
  n_samples <- nrow(newx)
  n_classes <- length(unique(unlist(strsplit(names(ovo_classifier_list), "_"))))
  
  # Initialize an empty matrix to store the probabilities
  e_hat <- matrix(0, nrow = n_samples, ncol = n_classes)
  colnames(e_hat) <- unique(unlist(strsplit(names(ovo_classifier_list), "_")))
  
  for (i in 1:n_classifiers) {
    # Predict probabilities using the i-th binary classifier
    # @Maren: change to ml_methods[[ml]]
    if (classifier == "svm") {
      probabilities <- predict(ovo_classifier_list[[i]], newdata = newx, probability = TRUE) %>% attr("probabilities")
    }
    # Compute probabilities for second class if not provided by default
    # @Maren: check how to add correct colnames
    if (ncol(probabilities==1)){
      probabilities[,2] <- 1 - probabilities[,1]
    }

    # Extract class names from the binary classifier name
    class_names <- unlist(strsplit(names(ovo_classifier_list)[i], "_"))
    class_i <- class_names[1]
    class_j <- class_names[2]
    
    # Update the corresponding columns in e_hat
    e_hat[, class_i] <- e_hat[, class_i] + probabilities[, class_i]
    e_hat[, class_j] <- e_hat[, class_j] + probabilities[, class_j]
  }
  
  # Normalize the probabilities to sum up to 1 for each sample
  e_hat <- e_hat / rowSums(e_hat)
  
  return(e_hat)
}


# Example usage:
# Generate some sample data
set.seed(42)

labels <- as.factor(W_training$W)

# Create the one-vs-one classifier
ovo_classifiers <- fit_ovo_classifier(X_training, labels, classifier = "svm")

# Access and use the binary classifiers as needed
e_hat_ovo <- predict_ovo_classifier(ovo_classifiers, newx = X_test)
bs_ovo <- brier_score(e_hat_ovo, W_test, binary = FALSE)
