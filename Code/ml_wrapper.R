#' Calculates arithmetic mean.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns list containing mean and number of observations
#'
#' @keywords internal
#'
mean_fit = function(x,y) {
  mean = mean(y)
  output = list("mean"=mean,"n"=nrow(x))
  output
}

#' Predicts arithmetic mean and provides prediction weights if required.
#' @param mean_fit Output of \code{\link{mean_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.mean_fit = function(mean_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) fit = rep(mean_fit$mean,nrow(x))
  else fit = rep(mean_fit$mean,nrow(xnew))
  
  if (isTRUE(weights)) w = matrix(1 / length(y),nrow(xnew),nrow(x))
  else w = NULL
  
  list("prediction"=fit, "weights"=w)
}


#' Calculates OLS fit.
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#'
#' @return Returns OLS coefficients
#'
#' @keywords internal
#'
ols_fit = function(x,y) {
  x = cbind(rep(1,nrow(x)),x)
  ols_coef = lm.fit(x,y)$coefficients
  ols_coef
}


#' Prediction based on OLS and provides prediction weights if required.
#' @param ols_fit Output of \code{\link{ols_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.ols_fit = function(ols_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  x = cbind(rep(1,nrow(x)),x)
  xnew = cbind(rep(1,nrow(xnew)),xnew)
  
  # Remove variables that were dropped due to collinearity
  x = x[,!is.na(ols_fit)]
  xnew = xnew[,!is.na(ols_fit)]
  
  # Calculate hat matrix
  hat_mat = xnew %*% solve(crossprod(x),tol=2.225074e-308) %*% t(x)
  fit = hat_mat %*% y
  
  if (weights==FALSE) hat_mat = NULL
  
  list("prediction"=fit,"weights"=hat_mat)
}


#' This function estimates cross-validated ridge regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
#' @keywords internal
#'
ridge_fit = function(x,y,args=list()) {
  ridge = do.call(cv.glmnet,c(list(x=x,y=y,alpha=0),args))
  ridge
}

#' Prediction based on Ridge and provides prediction weights if required.
#' @param ridge_fit Output of \code{\link{ridge_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.ridge_fit = function(ridge_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  fit = predict(ridge_fit,newx=xnew,type="response")
  
  if (weights==FALSE) hat_mat = NULL
  else {
    # Get covariate matrices
    n = nrow(x)
    p = ncol(x)
    x = scale(x)
    x = cbind(rep(1,nrow(x)),x)
    xnew = scale(xnew)
    xnew = cbind(rep(1,nrow(xnew)),xnew)
    
    # Calculate hat matrix, see also (https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat)
    hat_mat = xnew %*% solve(crossprod(x) + ridge_fit$lambda.min  * n / sd(y) * diag(x = c(0, rep(1,p)))) %*% t(x)
    fit = hat_mat %*% y
  }
  
  list("prediction"=fit,"weights"=hat_mat)
}


#' This function estimates cross-validated Post-Lasso based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "plasso"
#'
#' @keywords internal
#'
plasso_fit = function(x,y,args=list()) {
  plasso = do.call(plasso,c(list(x=x,y=y),args))
  plasso
}


#' Prediction based on Post-Lasso and provides prediction weights if required.
#' @param plasso_fit Output of \code{\link{plasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.plasso_fit = function(plasso_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  x = add_intercept(x)
  xnew = add_intercept(xnew)
  
  # Fitted values for post lasso
  nm_act = names(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl])[which(coef(plasso_fit$lasso_full)[,plasso_fit$ind_min_pl] != 0)]
  
  xact = x[,nm_act,drop=F]
  xactnew = xnew[,nm_act,drop=F]
  
  # Remove potentially collinear variables
  coef = lm.fit(xact,y)$coefficients
  xact = xact[,!is.na(coef)]
  xactnew = xactnew[,!is.na(coef)]
  
  hat_mat = xactnew %*% solve(crossprod(xact),tol=2.225074e-308) %*% t(xact)
  fit_plasso = hat_mat %*% y
  if (weights==FALSE) hat_mat = NULL
  
  list("prediction"=fit_plasso,"weights"=hat_mat)
}


#' Calculates Random Forest fit using the \code{\link{grf}} package
#'
#' @param x Matrix of covariates
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{regression_forest}}
#' @import grf
#'
#' @return An object with S3 class "regression_forest"
#'
#' @keywords internal
#'
forest_grf_fit = function(x,y,args=list()) {
  rf = do.call(regression_forest,c(list(X=x,Y=y),args))
  rf
}


#' Prediction based on Random Forest and provides prediction weights if required.
#' @param forest_grf_fit Output of \code{\link{regression_forest}} or \code{\link{forest_grf_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights If TRUE, weights underlying the prediction for xnew calculated
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{If \code{weights=TRUE} prediction weights of dimension \code{nrow(xnew)} x \code{nrow(x)}
#' containing the weights that deliver predictions where each row gives the weight that each training
#' outcome received in the prediction for xnew.}
#'
#' @keywords internal
#'
predict.forest_grf_fit = function(forest_grf_fit,x,y,xnew=NULL,weights=FALSE) {
  if (is.null(xnew)) xnew = x
  
  fit = predict(forest_grf_fit,newdata=xnew)$prediction
  
  if (weights==TRUE) {
    if (packageVersion("grf") < "2.0.0") w = get_sample_weights(forest_grf_fit,newdata=xnew)
    else  w = get_forest_weights(forest_grf_fit,newdata=xnew)
  }
  else w = NULL
  
  list("prediction"=fit, "weights"=w)
}




#' This function estimates cross-validated lasso regression based on the \code{\link{glmnet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes
#' @param args List of arguments passed to  \code{\link{glmnet}}
#' @import glmnet
#'
#' @return An object with S3 class "glmnet"
#'
#' @keywords internal
#'
lasso_fit = function(x,y,args=list()) {
  lasso = do.call(cv.glmnet,c(list(x=x,y=y),args))
  lasso
}


#' Prediction based on Lasso Forest and provides prediction weights if required.
#' @param lasso_fit Output of \code{\link{glmnet}} or \code{\link{lasso_fit}}
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param xnew Covariate matrix of test sample
#' @param weights Always FALSE as
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not available for Lasso, only for Post-Lasso}
#'
#' @keywords internal
#'
predict.lasso_fit = function(lasso_fit,x,y,xnew=NULL,weights=FALSE) {
  f = function() stop("No weighted representation of Lasso available.",call.=FALSE)
  if (isTRUE(weights)) f()
  if (is.null(xnew)) xnew = x
  
  fit = predict(lasso_fit,newx=xnew,type="response",s="lambda.min")
  
  list("prediction"=fit,"weights"="No weighted representation of Lasso available.")
}


################################################################################
################################################################################
############ Hier fangen meine Funktionen an ###################################
################################################################################
################################################################################



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

#' This function calculates the Bernulli Naive Bayes model based on the \code{\link{naivebayes}} package 
#'
#' @param x Covariate matrix of training sample
#' @param y Vector of outcomes of training sample
#' @param args List of arguments passed to  \code{\link{naivebayes}}
#' @import naivebayes
#'
#' @return An object with S3 method for class 'bernulli_naive_bayes'
#'
nb_bernulli_fit = function(x,y,args=list()){
  y = as.factor(y)
  model = do.call(naive_bayes, c(list(x=x,y=y),args))
  model
}
#' Prediction based on Bernulli Naive Bayes model.
#' @param nb_bernulli_fit Output of \code{\link{naive_bayes}} or \code{\link{nb_bernulli_fit}}
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
predict.nb_bernulli_fit = function(nb_bernulli_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  fit = predict(nb_bernulli_fit, newdata=xnew, type = "prob")
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
#xgboost_fit = function(x,y,args=list()){
#  K = length(unique(y))
#  if (K==2){
#     model = do.call(xgboost, c(list(data=x,label=y, nrounds=10, 
#                                     objective = "binary:logistic", eval_metric = "logloss", max_depth = 5),args)
#  }else{
#    y = y-1
#    model = do.call(xgboost, c(list(data=x,label=y, nrounds=10, num_class = K, objective = "multi:softprob", eval_metric = "mlogloss", max_depth = 5),args)
#  }
#  
#  model
#}

xgboost_fit <- function(x, y, args=list()) {
  K <- length(unique(y))
  
  if (K == 2) {
    model <- xgboost(data = x, label = y, nrounds = 10,
                     objective = "binary:logistic", eval_metric = "logloss", max_depth = 5, verbose=0, args)
  } else {
    if (min(y)==1) y <- y - 1
    model <- xgboost(data = x, label = y, nrounds = 10, num_class = K,
                     objective = "multi:softprob", eval_metric = "mlogloss", max_depth = 5, verbose=0, args)
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
#' @param args List of arguments passed to  \code{\link{xgboost}}
#' @import xgboost
#'
#' @return An object with S3 method for class 'svm'
#'

svm_fit = function(x,y,args=list()){
  model = do.call(svm, c(list(y = y, x = x, probability = TRUE, kernel = "sigmoid", 
                              type = 'C-classification', gamma = 0.1), args))
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
bagging_fit = function(x, y, args = list()) {
  data <- data.frame(y = as.factor(y), x)
  model = do.call(bagging, c(list(formula = y ~ ., data = data,
                                  mfinal = 50,
                                  control = rpart.control(objective = "binary:logistic",
                                                          eval_metric = "rmse",
                                                          eta = 0.3, max_depth = 2,
                                                          minsplit = 10,  # Adjust minsplit to an appropriate value
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


# #' Calculates AdaBoost.M1 fit using the \code{\link{adabag}} package
# #'
# #' @param x Matrix of covariates
# #' @param y vector of outcomes
# #' @param args List of arguments passed to  \code{\link{boosting}}
# #' @import adabag
# #'
# #' @return An object with S3 class "boosting"
# #'
# #' @keywords internal
# #'
# boosting_fit = function(x, y, args = list()) {
#   data <- data.frame(y = as.factor(y), x)
#   model = do.call(boosting, c(list(formula = y ~ ., data = data,
#                                   mfinal = 50,
#                                   control = rpart.control(objective = "binary:logistic",
#                                                           eval_metric = "rmse",
#                                                           eta = 0.3, max_depth = 2,
#                                                           minsplit = 10,  # Adjust minsplit to an appropriate value
#                                                           cp = -1
#                                   )), args))
#   model
# }
# #' Prediction based on AdaBoost.M1.
# #' @param boosting_fit Output of \code{\link{boosting}} or \code{\link{boosting_fit}}
# #' @param x Covariate matrix of training sample
# #' @param y Vector of outcomes of training sample
# #' @param xnew Covariate matrix of test sample
# #' @param weights Always FALSE
# #'
# #' @return Returns list containing:
# #' \item{prediction}{vector of predictions for xnew}
# #' \item{weights}{Not supported for propensity score estimation}
# #'
# #' @keywords internal
# #'
# predict.boosting_fit = function(boosting_fit,x,y,xnew=NULL,weights=FALSE){
#   if (is.null(xnew)) xnew = x
#   if (weights==TRUE) {
#     warning("Weights are not supported for propensity score estimation.")
#   }
#   
#   data = as.data.frame(xnew)
#   data$y <- as.factor(0)
#   
#   fit = predict(boosting_fit, newdata = data)$prob
#   if (length(unique(y))==2){
#     colnames(fit) = sort(unique(y))
#   }  
#   list("prediction"=fit,  "weights"="No weighted representation available.")
# }



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
adaboost_fit = function(x,y,args=list()){
  data <- data.frame(y = as.factor(y), x)
  model = do.call(adaboost ,c(list(data=data, formula = y ~ ., nIter=20, loss="logistic", w=NULL), args))
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
knn_fit = function(x,y,args=list()){ # args=list(kmax=floor(0.05*nrow(x)))
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
ranger_fit = function(x,y,args=list()){ #args=list(min.node.size=0.1*nrow(x), mtry=ceiling(sqrt(ncol(x))))
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
multinom_fit = function(x,y,args=list()){
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
#' \item{prediction}{vector of predictions for xnew}
#' \item{weights}{Not supported for propensity score estimation}
#'
#' @keywords internal
#'
predict.multinom_fit = function(multinom_fit, x,y,xnew=NULL,weights=FALSE){
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
ovo_fit <- function(x, y, classifier = "logit") {
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
      
      if (classifier=="xgboost"){
        subset_y <- ifelse(subset_y == i, 1, 0)
        
      }
      # @Maren: change to ml_methods[[ml]]
      ovo_classifiers[[paste(class_i, class_j, sep = "_")]] <- do.call(
        paste0(classifier, "_fit"),
        list(y = subset_y,
             x = subset_x)
      )
    }
  }
  
  return(ovo_classifiers)
}


#predict.ovo_fit <- function(ovo_fit,x,y, xnew=NULL,weights=FALSE, classifier = "logit") {
#  if (is.null(xnew)) xnew = x
#  if (weights==TRUE) {
#    warning("Weights are not supported for propensity score estimation.")
#  }
#  
#  n_classifiers <- length(ovo_fit)
#  n_samples <- nrow(xnew)
#  n_classes <- length(unique(unlist(strsplit(names(ovo_fit), "_"))))
#  
#  # Initialize an empty matrix to store the probabilities
#  fit <- matrix(0, nrow = n_samples, ncol = n_classes) %>% as_tibble()
#  colnames(fit) <- unique(unlist(strsplit(names(ovo_fit), "_")))
#  
#  for (i in 1:n_classifiers) {
#    # Extract class names from the binary classifier name
#    class_names <- unlist(strsplit(names(ovo_fit)[[i]], "_"))
#    class_i <- class_names[1]
#    class_j <- class_names[2]
#    subset_y <- y[y %in% c(class_i, class_j)]
#    #print(paste0("classifier ", i))
# 
#    # Predict probabilities using the i-th binary classifier
#    fit_raw <- do.call(paste0("predict.", classifier,  "_fit"), 
#                       list(ovo_fit[[i]], x=x, y=subset_y, xnew = xnew))
#
#
#    # Compute probabilities for second class if not provided by default
#    # @Maren: check how to add correct colnames
#    fit_raw = fit_raw %>% as_tibble()
#    if (ncol(fit_raw==1)){
#      fit_raw[,2] <- 1 - fit_raw[,1]
#    }
#    if (!(all(class_names %in% colnames(fit_raw)))) {
#      colnames(fit_raw) <- class_names
#    }
#    
#    
#    # Update the corresponding columns in e_hat
#    fit[, `class_i`] <- fit[, `class_i`] + fit_raw[, `class_i`]
#    fit[, `class_j`] <- fit[, `class_j`] + fit_raw[, `class_j`]
#  }
#  
#  # Normalize the probabilities to sum up to 1 for each sample
#  ##@Maren: research
#  fit <- fit / rowSums(fit)
#  
#  
#  #list("prediction"=fit,  "weights"="No weighted representation available.")
#  fit
#}


predict.ovo_fit <- function(ovo_fit,x,y, xnew=NULL,weights=FALSE, classifier = "logit") {
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovo_fit)
  n_samples <- nrow(xnew)
  n_classes <- length(unique(unlist(strsplit(names(ovo_fit), "_"))))
  
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
    fit_raw <- do.call(paste0("predict.", classifier,  "_fit"), 
                       list(ovo_fit[[i]], x=x, y=subset_y, xnew = xnew))
    
    
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
  
  
  # Assign the results to the fit matrix
  fit <- t(simplify2array(opt_results))
  colnames(fit) <- unique(unlist(strsplit(names(ovo_fit), "_")))
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}



# One-vs-Rest Classifier
ovr_fit <- function(x, y, classifier = "logit") {
  class_labels <- sort(unique(y))
  num_classes <- length(class_labels)
  ovr_classifiers <- list()
  
  # Create binary classifiers for each pair of classes
  for (i in 1:(num_classes)) {
      class_i <- class_labels[i]      
      binarized_y <- ifelse(y == i, 1, 0)
  
      
      # @Maren: change to ml_methods[[ml]]
      ovr_classifiers[[paste(class_i)]] <- do.call(
        paste0(classifier, "_fit"),
        list(y = binarized_y, x = x)
      )
  }
  
  return(ovr_classifiers)
}


predict.ovr_fit <- function(ovr_fit,x,y, xnew=NULL,weights=FALSE, classifier = "logit") {
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    warning("Weights are not supported for propensity score estimation.")
  }
  
  n_classifiers <- length(ovr_fit)
  n_samples <- nrow(xnew)
  n_classes <- n_classifiers

  # Initialize an empty matrix to store the probabilities
  fit <- matrix(0, nrow = n_samples, ncol = n_classes) %>% as_tibble()
  colnames(fit) <- seq.int(1, n_classes)
  
  for (i in 1:n_classifiers) {
    binarized_y <- ifelse(y == i, 1, 0)
    #print(paste0("classifier ", i))
    # Predict probabilities using the i-th binary classifier
    # @Maren: change to ml_methods[[ml]]
    fit_raw <- do.call(paste0("predict.", classifier,  "_fit"), 
                       list(model[[i]], x=x, y=binarized_y, xnew = xnew))$prediction

    fit_raw = as_tibble(fit_raw) %>% select("1") 


    # Update the corresponding columns in e_hat
    fit[ , `i` ] <- fit[ , `i` ] + fit_raw
  }
  
  # Normalize the probabilities to sum up to 1 for each sample
  fit <- fit / rowSums(fit)
  
  
  #list("prediction"=fit,  "weights"="No weighted representation available.")
  fit
}










#model_fit = function(x,y,args=list()){
#  model = do.call( ... ,c(list(...), args))
#  model
#}
#
#predict.model_fit = function(model_fit, x,y,xnew=NULL,weights=FALSE){
#  if (is.null(xnew)) xnew = x
#  if (weights==TRUE) {
#    warning("Weights are not supported for propensity score estimation.")
#  }
#  
#  fit = predict(model_fit, )
#  list("prediction"=fit,  "weights"="No weighted representation available.")
#}