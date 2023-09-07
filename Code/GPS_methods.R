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
source("MetaLearners_tools.R");rm(cl, nc)
set.seed(42)


# create toy dataset
outcome <- function(t, X){
  return( (1+t)*X )
}

n <- 5000
K <- 5

X1 <- runif(n = n)
X2 <- runif(n = n)
X3 <- runif(n = n)
W <- sample(seq(0, 1, by = 1/(K-1)), replace = TRUE, n)

W.levels <- sort(unique(W))
K <- length(W.levels)

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
r_hat_unif <- matrix(NA, nrow = nrow(test_data), ncol = K)
for(k in 1:K){
  r_hat_unif[,k] = apply(data.frame(W.levels[k], test_data$X1), 1, K_GPS_unif)
}

# 1. logistic regression
r_hat_logit <- matrix(NA, nrow = nrow(test_data), ncol = K)
# model_logit = mlogit::mlogit(W ~ X, data = training_data) # need the data to be reshaped 
model_logit <- nnet::multinom(formula, data = training_data, model = TRUE)
r_hat_logit <- predict(model_logit, newdata = X_test, type = "prob")

# 2. logistic regression with cross-validation
r_hat_logit_cv <- matrix(NA, nrow = nrow(test_data), ncol = K)
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
#r_hat_logit_cv <- predict(model_logit_cv, newdata = X_test, type = "prob")

# 3. Ridge classification
model_ridge <- glmnet(x = as.matrix(X_training), y = as.matrix(W_training), alpha = 0, family = "multinomial")
r_hat_ridge <- predict(model_ridge, newx = as.matrix(X_test), s = 0.01, type = "response")

# 4. Ridge classification cross-validation
model_ridge_cv <- cv.glmnet(x = as.matrix(X_training), y = as.matrix(W_training), alpha = 0, family = "multinomial")
r_hat_ridge_cv <- predict(model_ridge_cv, newx = as.matrix(X_test), s = 0.01, type = "response")

# 5. Naive bayes classifier (Bernulli)
model_nb_bernulli <- bernoulli_naive_bayes(X_training, as.factor(training_data$W))
r_hat_nb_bernulli <- predict(model_nb_bernulli, newdata = as.matrix(X_test), type = "prob")

# 6. Naive bayes classifier (Gaussian)
model_nb_gaussian <- bernoulli_naive_bayes(X_training, as.factor(training_data$W))
r_hat_nb_gaussian <- predict(model_nb_gaussian, newdata = as.matrix(X_test), type = "prob")

all(r_hat_nb_bernulli == r_hat_nb_gaussian)

# 7. Decision Tree classifier
model_decision_tree <- rpart(formula, data = training_data, method = "class")
r_hat_decision_tree <- predict(model_decision_tree, X_test, type = "prob")

# 8. Extra Tree classifier
model_extra_tree <-
r_hat_extra_tree <-

# 9. Random Forest
r_hat_rf <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_rf <- randomForest(y = as.factor(training_data$W), x = X_training)
r_hat_rf <- predict(model_rf, X_test, type = "prob")

# 10. Label Propagation
model_label_prop <-
r_hat_label_prop <-

# 11. Label Spreading
model_label_spread <-
r_hat_label_spread <-

# 12. Linear Support Vector Machine
r_hat_svm <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_svm <- svm(y = as.factor(training_data$W), x = X_training, probability = TRUE)
r_hat_svm <- predict(model_svm, X_test, probability = TRUE) %>% attr("probabilities")
# 13. K-Nearest-Neighbors
r_hat_knn <- matrix(NA, nrow = nrow(test_data), ncol = K)
model_knn <- FNN::knn(train = as.matrix(X_training), test = as.matrix(X_test), cl = training_data$W, prob = TRUE)
model_knn <- class::knn(train = as.matrix(X_training), test = as.matrix(X_test), cl = training_data$W, prob = TRUE)
r_hat_knn <- attributes(model_knn)$prob

# 14. Nearest Centroid
model_nearest_centriod <-
r_hat_nearest_centriod <-

# 15. Nearest-Neighbor (Radius)
model_nn_radius <-
r_hat_nn_radius <-

# 16. Linear Discriminant Analysis
model_lda <- lda(formula, data = training_data)
r_hat_lda <- predict(model_lda, X_test, type = "prob")$posterior

# 17. Quadratic Discriminant Analysis
model_qda <- qda(formula, data = training_data)
r_hat_qda <- predict(model_qda, X_test, type = "prob")$posterior

# 18. Multi-layer Perceptron classifier
model_mlpc <- nnet(x = scale(X_training), y = class.ind(training_data$W), linout = FALSE, size = c(10, 5), softmax = TRUE, maxit = 200)
r_hat_mlpc <- predict(model_mlpc, X_test)

# 19. Stochastic Gradient Descent (SGD)
model_sgd <-
r_hat_sgd <-
  