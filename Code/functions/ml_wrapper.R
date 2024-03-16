#' This function calculates the (multinomial) logistic regression based on the \code{\link{glmnet}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
logit_fit = function(x,y,args=list()){
  if (length(unique(y))==2){
    logit = do.call(glmnet, c(list(x=x,y=y, family = "binomial", type.measure = "class"),args))
  }else{
    logit = do.call(glmnet, c(list(x=x,y=y, family = "multinomial", type.measure = "class"),args))
  }
  logit
}


#' Prediction based on (multinomial) logistic regression.
#' @param logit_fit Output of \code{\link{glmnet}} or \code{\link{logit_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.logit_fit = function(logit_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  lambda = min(logit_fit$lambda)
  fit = predict(logit_fit, newx=as.matrix(xnew), s = lambda, type = "response") %>% as_tibble()
  if (length(logit_fit$glmnet.fit$classnames)==2){
    fit[,2] = 1 - fit[,1]
    colnames(fit) = rev(logit_fit$glmnet.fit$classnames)
      }else{
    colnames(fit) = logit_fit$classnames
  }
  
  fit
}


#' This function calculates the Gaussian Naive Bayes model based on the \code{\link{naivebayes}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{naivebayes}}
#' @import naivebayes
#'
#' @return An object with S3 method for class 'gaussian_naive_bayes'
#'
nb_gaussian_fit = function(x,y,args=list()){
  y = as.factor(y)
  model = do.call(gaussian_naive_bayes, c(list(x=x,y=y),args))
  model
}

#' Prediction based on Gaussian Naive Bayes model.
#' @param nb_gaussian_fit Output of \code{\link{gaussian_naive_bayes}} or \code{\link{nb_gaussian_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.nb_gaussian_fit = function(nb_gaussian_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(nb_gaussian_fit, newdata=xnew, type = "prob")
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}

#' This function calculates the bernoulli Naive Bayes model based on the \code{\link{naivebayes}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{naivebayes}}
#' @import naivebayes
#'
#' @return An object with S3 method for class 'bernoulli_naive_bayes'
#'
nb_bernoulli_fit = function(x,y,args=list()){
  y = as.factor(y)
  model = do.call(naive_bayes, c(list(x=x,y=y),args))
  model
}
#' Prediction based on bernoulli Naive Bayes model.
#' @param nb_bernoulli_fit Output of \code{\link{naive_bayes}} or \code{\link{nb_bernoulli_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.nb_bernoulli_fit = function(nb_bernoulli_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(nb_bernoulli_fit, newdata=xnew, type = "prob")
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}



#' This function calculates the xgboost model based on the \code{\link{xgboost}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{xgboost}}
#' @import xgboost
#'
#' @return An object with S3 method for class 'xgb.Booster'
#'

# seems not to work with the do.call() function
xgboost_fit = function(x,y,args=list(nrounds=40, eta=0.01)){
  K = length(unique(y))
  if (K==2){
     model = do.call(xgboost, c(list(data=x,label=y, verbose=0,
                                     objective = "binary:logistic", eval_metric = "logloss"),args))
  }else{
    if (min(y)==1) y <- y - 1
    model = do.call(xgboost, c(list(data=x,label=y, num_class = K, verbose=0,
                                    objective = "multi:softprob", eval_metric = "mlogloss"),args))
  }
  
  model
}



#' Prediction based on xgboost model.
#' @param xgboost_fit Output of \code{\link{xgboost}} or \code{\link{xgboost_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.xgboost_fit = function(xgboost_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(xgboost_fit, newdata=xnew, type = "prob") 
  if (length(unique(y))==2){
  #if (length((unique(subset_y)))==2){
    fit = as_tibble(fit)
    fit[,2] = 1 - fit[,1]
    colnames(fit) = rev(sort(unique(y)))
    #colnames(fit) = rev(sort(unique(subset_y)))
    
  }else{
    fit = matrix(fit, nrow = nrow(xnew), ncol = length(unique(y)), byrow = TRUE)
    colnames(fit) = sort(unique(y))
    
  }

  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}


#' This function calculates the Support Vector Machines model based on the \code{\link{e1071}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{svm}}
#' @import e1071
#'
#' @return An object with S3 method for class 'svm'
#'

svm_fit = function(x,y,args=list(gamma = 0.1)){
  model = do.call(svm, c(list(y = y, x = x, probability = TRUE, kernel = "radial", 
                              type = 'C-classification'), args))
  return(model)
}

#' Prediction based on Support Vector Machines model.
#' @param svm_fit Output of \code{\link{e1071::svm}} or \code{\link{svm_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.svm_fit <- function(svm_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  # Perform predictions using the trained SVM model
  fit <- predict(svm_fit, newdata = xnew, probability = TRUE) %>% attr("probabilities")
  
  # Return the predictions and weights (if any)
  result <- list("prediction" = fit, "weights"="No weighted representation available.")
  return(result)
}



#' This function calculates the Linear Discriminant Analysis based on the \code{\link{MASS}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{MASS}}
#' @import MASS
#'
#' @return An object with S3 method for class 'lda'
#'
lda_fit = function(x,y,args=list()){
  model = do.call(lda, c(list(x = x, grouping = y, cv = TRUE)))
  model
}

#' Prediction based on Linear Discriminant Analysis.
#' @param lda_fit Output of \code{\link{lda}} or \code{\link{lda_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.lda_fit = function(lda_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  fit <- predict(lda_fit, xnew, type="prob")$posterior # @Maren: remove type and retry :)
  
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}

#' This function calculates the Quadratic Discriminant Analysis based on the \code{\link{MASS}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{MASS}}
#' @import MASS
#'
#' @return An object with S3 method for class 'qda'
#'
qda_fit = function(x,y,args=list()){
  y = as.factor(y)
  
  model = do.call(qda, c(list(x = x, grouping = y, cv = TRUE), args))
  model
}

#' Prediction based on Quadratic Discriminant Analysis.
#' @param qda_fit Output of \code{\link{qda}} or \code{\link{qda_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.qda_fit = function(qda_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit <- predict(qda_fit, xnew, type="prob")$posterior
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}


#' Calculates Probability Forest fit using the \code{\link{grf}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{probability_forest}}
#' @import grf
#'
#' @return An object with S3 class "probability_forest"
#'
#' @keywords internal
#'
probability_forest_fit = function(x,y,args=list(num.trees=1000)){
  model = do.call(probability_forest, c(list(X = x, Y = as.factor(y)), args))
  model
}
#' Prediction based on Probability Forest fit.
#' @param probability_forest_fit Output of \code{\link{probability_forest}} or \code{\link{probability_forest_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.probability_forest_fit = function(probability_forest_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  
  if (weights==TRUE) {
    if (packageVersion("grf") < "2.0.0") w = get_sample_weights(forest_grf_fit,newdata=xnew)
    else  w = get_forest_weights(forest_grf_fit,newdata=xnew)
  }
  else w = NULL
  
  fit = predict(probability_forest_fit, newdata=xnew)$predictions
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}


#' Calculates Bagged classification tree model using the \code{\link{adabag}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{bagging}}
#' @import adabag
#'
#' @return An object with S3 class "bagging"
#'
#' @keywords internal
#'
bagging_fit = function(x, y, args = list(eta = 0.3, 
                                         max_depth = 2,
                                         minsplit = 10)) {
  data <- data.frame(y = as.factor(y), x)
  model = do.call(bagging, c(list(formula = y ~ ., data = data,
                                  mfinal = 50,
                                  control = rpart.control(objective = "binary:logistic",
                                                          eval_metric = "rmse",
                                                          cp = -1
                                  )), args))
  model
}

#' Prediction based on Bagged classification tree model.
#' @param bagging_fit Output of \code{\link{bagging}} or \code{\link{bagging_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.bagging_fit = function(bagging_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  data = as.data.frame(xnew)
  data$y <- as.factor(0)
  
  fit = predict(bagging_fit, newdata = data)$prob
  if (length(unique(y))==2){
    colnames(fit) = sort(unique(y))
  }
  list("prediction"=fit,  "weights"="No weighted representation available.")
}



#' Calculates Adaboost.M1 model using the \code{\link{fastadaboost}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{adaboost}}
#' @import fastAdaboost
#'
#' @return An object S3 method for class 'adaboost'
#'
#' @keywords internal
#'
adaboost_fit = function(x,y,args=list(nIter=20)){
  data <- data.frame(y = as.factor(y), x)
  model = do.call(adaboost ,c(list(data=data, formula = y ~ ., loss="logistic", w=NULL), args))
  model
}
#' Prediction based on Adaboost.M1 model.
#' @param adaboost_fit Output of \code{\link{adaboost}} or \code{\link{adaboost_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.adaboost_fit = function(adaboost_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(adaboost_fit, newdata=xnew, type = "probs")$prob
  if (length(unique(y))==2){
    colnames(fit) = sort(unique(y))
  }
  list("prediction"=fit,  "weights"="No weighted representation available.")
}


#' Calculates k-Nearest Neighbor model using the \code{\link{kknn}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{kknn}}
#' @import grf
#'
#' @return An object with S3 class "kknn"
#'
#' @keywords internal
#'
# @Maren: kmax depedent on # classes?
knn_fit = function(x,y,args=list(kmax=floor(0.05*nrow(x)))){ 
  data <- data.frame(y = as.factor(y), x)
  model = do.call(train.kknn ,c(list(formula = y ~ ., data = data), args)) # args
  model
}

#' Prediction based on k-Nearest Neighbor model.
#' @param knn_fit Output of \code{\link{kknn}} or \code{\link{knn_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.knn_fit = function(knn_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  xnew = as.data.frame(xnew)
  fit = predict(knn_fit, newdata = xnew, type='prob')
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}



#' Calculates k-Nearest Neighbor model using the \code{\link{kknn}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{kknn}}
#' @import grf
#'
#' @return An object with S3 class "kknn"
#'
#' @keywords internal
#'
knn_radius_fit = function(x,y,args=list(distance=10)){
  data <- data.frame(y = as.factor(y), x)
  model = do.call(train.kknn ,c(list(formula = y ~ ., data = data), args))
  model
}
#' Prediction based on k-Nearest Neighbor model.
#' @param knn_radius_fit Output of \code{\link{kknn}} or \code{\link{knn_radius_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.knn_radius_fit = function(knn_radius_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  xnew = as.data.frame(xnew)
  fit = predict(knn_radius_fit, newdata = xnew, type='prob')
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}



#' Calculates Multilayer Perceptron model using the \code{\link{nnet}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{nnet}}
#' @import nnet
#'
#' @return An object with S3 class "nnet"
#'
#' @keywords internal
#'
mlpc_fit = function(x,y,args=list(size = c(5))){ #size=1
  #if(is.null(size)) size=1
  data <- data.frame(y = as.factor(y), x)
  model = do.call(nnet ,c(list(x=x, y = class.ind(y), softmax = TRUE, lineout=TRUE), args))
  #model = do.call(neuralnet, c(list(formula = y ~ ., data = data), args))
  model
}
#' Prediction based on Multilayer Perceptron fit.
#' @param nnet_fit Output of \code{\link{nnet}} or \code{\link{nnet_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.mlpc_fit = function(mlpc_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(mlpc_fit, newdata=xnew)
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}


#' Calculates BART classification model using the \code{\link{bartMachine}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{bartMachine}}
#' @import bartMachine
#'
#' @return An object with S3 class "bartMachine"
#'
#' @keywords internal
#'
bart_fit = function(x,y,args=list()){
  model = do.call(bartMachine, c(list(X=as.data.frame(x), y=as.factor(y)), args))
  model
}
#' Prediction based on Probability Forest fit.
#' @param bart_fit Output of \code{\link{bartMachine}} or \code{\link{bart_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.bart_fit = function(bart_fit, x,y,xnew, weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(bart_fit, new_data=as.data.frame(xnew)) 
  fit = as_tibble(fit)
  fit[,2] = 1 - fit[,1]
  colnames(fit) = rev(sort(unique(y)))
  list("prediction"=fit,  "weights"="No weighted representation available.")
}

#' Calculates Probability Forest fit using the \code{\link{ranger}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{ranger}}
#' @import ranger
#'
#' @return An object with S3 class "ranger"
#'
#' @keywords internal
#'
ranger_fit = function(x,y,args=list(num.trees = 1000,
                                    min.node.size=0.1*nrow(x), 
                                    mtry=ceiling(sqrt(ncol(x))))){ #args=list()
  data <- data.frame(y = as.factor(y), x)

  model = do.call(ranger, c(list(data=data, formula= y~., probability=TRUE), args))
  model
}

#' Prediction based on the \code{\link{ranger}} Probability Forest fit.
#' @param probability_forest_fit Output of \code{\link{ranger}} or \code{\link{probability_forest_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ranger_fit = function(ranger_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  data = as.data.frame(xnew)
  data$y <- as.factor(0)
  
  fit <- predict(ranger_fit, data = data)$predictions 

  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}

#model = ranger_fit(x = as.matrix(X_training), y = W_training)
#preds = predict.ranger_fit(model, x = as.matrix(X_training), y = W_training, xnew = X_test)$prediction

#' This function calculates the (multinomial) logistic regression based on the \code{\link{nnet}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{nnet}}
#' @import nnet
#'
#' @return An object with S3 method for class 'nnet'
#'
logit_nnet_fit = function(x,y,args=list()){
  data <- data.frame(y = as.factor(y), x)
  
  model = do.call(multinom, c(list(data=data, formula= y~., model=TRUE), args))
  model
}
#' Prediction based on (multinomial) logistic regression.
#' @param multinom_fit Output of \code{\link{nnet}} or \code{\link{multinom_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.logit_nnet_fit = function(multinom_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  data = as.data.frame(xnew)
  data$y <- as.factor(0)
  
  fit <- predict(multinom_fit, newdata = xnew,type="probs") %>% as_tibble()
  if (length(unique(y))==2){
    fit[,2] = 1-fit[,1]
    colnames(fit) = rev(sort(unique(y)))
  }
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}






# Function to create a one-vs-one classifier
#' Fits One-Vs-One (OvO) binary classifiers for multi-class classification
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost")
#' 
#' @return List of fitted OvO binary classifiers
#'
ovo_fit <- function(x, y, method = "logit") {
  
  if (min(y)==0) y = y+1
  
  class_labels <- sort(unique(y))
  num_classes <- length(class_labels)
  ovo_classifiers <- list()
  
  
  
  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes - 1)) {
    for (j in (i + 1):num_classes) {
      class_i <- class_labels[i]
      class_j <- class_labels[j]
      
      subset_x <- x[y %in% c(class_i, class_j),]
      subset_y <- y[y %in% c(class_i, class_j)]
      
      if (method=="xgboost"){
        subset_y <- ifelse(subset_y == i, 1, 0)
        
      }
      # @Maren: change to ml_methods[[ml]]
      ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
        paste0(method, "_fit"),
        list(y = subset_y,
             x = subset_x)
      )
    }
  }
  
  return(ovo_classifiers)
}


# Register parallel backend
# Adjust the number of cores as per your system configuration

# Function to create a one-vs-one classifier
#' Fits One-Vs-One (OvO) binary classifiers for multi-class classification in parallel
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost")
#' 
#' @return List of fitted OvO binary classifiers
#'
ovo_fit_parallel <- function(x, y, method = "logit") {
  
  # fist class must be class 1
  if (min(y) == 0) y <- y + 1
  
  class_labels <- sort(unique(y))
  num_classes <- length(class_labels)
  ovo_classifiers <- list()
  
  # Create binary classifiers for each pair of classes in parallel
  ovo_classifiers <- foreach(i = 1:(num_classes - 1), .combine = "c") %dopar% {
    classifiers <- list()
    for (j in (i + 1):num_classes) {
      class_i <- class_labels[i]
      class_j <- class_labels[j]
      
      subset_x <- x[y %in% c(class_i, class_j), ]
      subset_y <- y[y %in% c(class_i, class_j)]
      
      if (method == "xgboost") {
        subset_y <- ifelse(subset_y == i, 1, 0)
      }
      
      classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
        paste0(method, "_fit"),
        list(y = subset_y,
             x = subset_x)
      )
    }
    classifiers
  }
  
  return(ovo_classifiers)
}

# Now you can use ovo_fit_parallel instead of ovo_fit
#' Predicts class probabilities using One-Vs-One (OvO) fitted classifiers 
#' @param ovo_fit Output list of \code{\link{ovo_fit}} or \code{\link{ovo_fit_parallel}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Flag indicating whether weights are used (unsupported for propensity score estimation)
#' @param method Method used for fitting the binary classifiers (e.g., "logit" or "xgboost")
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ovo_fit <- function(ovo_fit,x,y, xnew=NULL,weights=FALSE, method = "logit") {
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovo_fit)
  n_samples <- nrow(xnew)
  n_classes <- length(unique(unlist(strsplit(names(ovo_fit), "_"))))
  
  if (min(y)==0) y = y+1
  
  # Create an empty q_matrix tensor
  q_matrix_tensor <- array(NA, dim = c(n_samples, n_classes, n_classes))
  
  # Populate q_matrix_tensor using vectorized operations
  for (i in 1:n_classifiers) {
    class_names <- unlist(strsplit(names(ovo_fit)[[i]], "_"))
    class_i <- class_names[1]
    class_j <- class_names[2]
    
    subset_y <- y[y %in% c(class_i, class_j)]

    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(paste0("predict.", method,  "_fit"), 
                       list(ovo_fit[[i]], x=x, y=subset_y, xnew = xnew))
    
    if(is.list(fit_raw)) if(!is.data.frame(fit_raw)) fit_raw = fit_raw$prediction
    # Compute probabilities for second class if not provided by default
    fit_raw = fit_raw %>% as_tibble()
    if (ncol(fit_raw==1)){
      fit_raw[,2] <- 1 - fit_raw[,1]
    }
    if (!(all(class_names %in% colnames(fit_raw)))) {
      colnames(fit_raw) <- class_names
    }
    
    q_matrix_tensor[, as.numeric(class_i),as.numeric(class_j)] <- fit_raw[, `class_i`][[1]] 
    q_matrix_tensor[, as.numeric(class_j),as.numeric(class_i)] <- fit_raw[, `class_j`][[1]] 
  }
  
  # Perform optimization for all samples simultaneously
  system.time({
    opt_results <- apply(q_matrix_tensor, MARGIN = 1, function(q_matrix_row) {
      opt_result <- optim(rep(1/n_classes, n_classes), 
                          kl_convergence, 
                          q_matrix = matrix(q_matrix_row, nrow = n_classes),
                          method = "L-BFGS-B", 
                          lower = rep(0, n_classes), upper = rep(1, n_classes))
      
      # Normalize p to sum to 1
      p_optimized <- opt_result$par / sum(opt_result$par)
      return(p_optimized)
    })
  })
  
  
  
  # Assign the results to the fit matrix
  fit <- t(simplify2array(opt_results))
  colnames(fit) <- unique(unlist(strsplit(names(ovo_fit), "_")))
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  return(fit)
}

# Now you can use ovo_fit_parallel instead of ovo_fit
#' Predicts class probabilities using One-Vs-One (OvO) fitted classifiers in parallel
#' @param ovo_fit Output list of \code{\link{ovo_fit}} or \code{\link{ovo_fit_parallel}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Flag indicating whether weights are used (unsupported for propensity score estimation)
#' @param method Method used for fitting the binary classifiers (e.g., "logit" or "xgboost")
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ovo_fit_parallel <- function(ovo_fit, x, y, xnew = NULL, weights = FALSE, method = "logit") {
  if (is.null(xnew)) xnew <- x
  if (weights) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovo_fit)
  n_samples <- nrow(xnew)
  n_classes <- length(unique(unlist(strsplit(names(ovo_fit), "_"))))
  
  if (min(y) == 0) y <- y + 1
  
  # Create an empty q_matrix tensor
  q_matrix_tensor <- array(NA, dim = c(n_samples, n_classes, n_classes))
  
  # Populate q_matrix_tensor using vectorized operations
  for (i in 1:n_classifiers) {
    class_names <- unlist(strsplit(names(ovo_fit)[[i]], "_"))
    class_i <- class_names[1]
    class_j <- class_names[2]
    
    subset_y <- y[y %in% c(class_i, class_j)]
    #print(paste0("classifier ", i))
    
    # Predict probabilities using the i-th binary classifier
    fit_raw <- do.call(paste0("predict.", method,  "_fit"), 
                       list(ovo_fit[[i]], x=x, y=subset_y, xnew = xnew))
    
    if(is.list(fit_raw)) if(!is.data.frame(fit_raw)) fit_raw = fit_raw$prediction
    # Compute probabilities for second class if not provided by default
    # @Maren: check how to add correct colnames
    fit_raw = fit_raw %>% as_tibble()
    if (ncol(fit_raw==1)){
      fit_raw[,2] <- 1 - fit_raw[,1]
    }
    if (!(all(class_names %in% colnames(fit_raw)))) {
      colnames(fit_raw) <- class_names
    }
    
    q_matrix_tensor[, as.numeric(class_i),as.numeric(class_j)] <- fit_raw[, `class_i`][[1]] 
    q_matrix_tensor[, as.numeric(class_j),as.numeric(class_i)] <- fit_raw[, `class_j`][[1]] 
  }
  
  # Perform optimization for all samples simultaneously
  system.time({
    opt_results <- foreach(q_matrix_row = 1:n_samples, .combine = "c") %dopar% {
      opt_result <- optim(rep(1/n_classes, n_classes), 
                          kl_convergence, 
                          q_matrix = q_matrix_tensor[q_matrix_row,,],
                          method = "L-BFGS-B", 
                          lower = rep(0, n_classes), upper = rep(1, n_classes))
      
      # Normalize p to sum to 1
      p_optimized <- opt_result$par / sum(opt_result$par)
      return(p_optimized)
    }
  })
  fit <- matrix(opt_results, nrow = n_samples, ncol = n_classes, byrow = TRUE)
  
  # Assign the results to the fit matrix
  colnames(fit) <- unique(unlist(strsplit(names(ovo_fit), "_")))
  
  return(fit)
}


# One-vs-Rest Classifier
#' Fits One-Vs-Rest (OvR) binary classifiers for multi-class classification 
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param method Method to use for fitting binary classifiers (e.g., "logit" or "xgboost")
#' 
#' @return List of fitted OvR binary classifiers
#'
ovr_fit <- function(x, y, method = "logit") {
  class_labels <- sort(unique(y))
  num_classes <- length(class_labels)
  ovr_classifiers <- list()
  
  if (min(y)==0) y = y+1
  
  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes)) {
      class_i <- class_labels[i]      
      binarized_y <- ifelse(y == i, 1, 0)
  
      
      # @Maren: change to ml_methods[[ml]]
      ovr_classifiers[[paste(class_i)]] <- do.call(
        paste0(method, "_fit"),
        list(y = binarized_y, x = x)
      )
  }
  
  return(ovr_classifiers)
}

# Now you can use ovo_fit_parallel instead of ovo_fit
#' Predicts class probabilities using One-Vs-One (OvR) fitted classifiers 
#' @param ovo_fit Output list of \code{\link{ovr_fit}} 
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Flag indicating whether weights are used (unsupported for propensity score estimation)
#' @param method Method used for fitting the binary classifiers (e.g., "logit" or "xgboost")
#'
#' @return Returns list containing:
#' \item{prediction}{matrix of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.ovr_fit <- function(ovr_fit,x,y, xnew=NULL,weights=FALSE, method = "logit") {
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovr_fit)
  n_samples <- nrow(xnew)
  n_classes <- n_classifiers
  
  if (min(y)==0) y = y+1
  
  # Initialize an empty matrix to store the probabilities
  fit <- matrix(0, nrow = n_samples, ncol = n_classes) %>% as_tibble()
  colnames(fit) <- seq.int(1, n_classes)
  
  for (i in 1:n_classifiers) {
    binarized_y <- ifelse(y == i, 1, 0)
    # Predict probabilities using the i-th binary classifier
    # @Maren: change to ml_methods[[ml]]
    fit_raw <- do.call(paste0("predict.", method,  "_fit"), 
                       list(ovr_fit[[i]], x=x, y=binarized_y, xnew = xnew))#$prediction
    #print(colnames(fit_raw))
    if(is.list(fit_raw)) if(!is.data.frame(fit_raw)) fit_raw = fit_raw$prediction
    if (ncol(fit_raw)==2){
      fit_raw = as_tibble(fit_raw) %>% select("1") 
    }else{
      fit_raw <- as_tibble(fit_raw)
    }


    # Update the corresponding columns in e_hat
    fit[ , `i` ] <- fit[ , `i` ] + fit_raw
  }
  
  # Normalize the probabilities to sum up to 1 for each sample
  fit <- fit / rowSums(fit)
  
  
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}




