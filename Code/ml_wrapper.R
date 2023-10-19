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
  
  list("prediction"=fit,"weights"=w)
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
  
  list("prediction"=fit,"weights"=w)
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


logit_fit = function(x,y,args=list()){
  logit = do.call(cv.glmnet, c(list(x=x,y=y),args))
  logit
}

predict.logit_fit = function(logit_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(logit_fit, newx=xnew, type = "response")
  list("prediction"=fit,"weights"=w)
}



model = logit_fit(as.matrix(X_training), as.matrix(W_training))
preds = predict.logit_fit(model, as.matrix(X_training), as.matrix(W_training), xnew=as.matrix(X_test))




nb_gaussian_fit = function(x,y,args=list()){
  y = as.factor(y)
  model = do.call(naive_bayes, c(list(x=x,y=y),args))
  model
}

predict.nb_gaussian_fit = function(nb_gaussian_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(nb_gaussian_fit, newdata=xnew, type = "prob")
  list("prediction"=fit,"weights"=w)
}

nb_bernulli_fit = function(x,y,args=list()){
  y = as.factor(y)
  model = do.call(bernoulli_naive_bayes, c(list(x=x,y=y),args))
  model
}

predict.nb_bernulli_fit = function(nb_bernulli_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(nb_bernulli_fit, newdata=xnew, type = "prob")
  list("prediction"=fit,"weights"=w)
}


model = nb_bernulli_fit(x=as.matrix(X_training), y=as.matrix(training_data$W))
preds = predict.nb_bernulli_fit(model, x = as.matrix(X_training), y = as.factor(training_data$W), xnew = as.matrix(X_test))$prediction





xgboost_fit = function(x,y,args=list()){
  model = do.call(xgboost, c(list(data=x,label=y, nrounds=10, objective = "binary:logistic", eval_metric = "rmse", max_depth = 5),args))
  model
}
predict.xgboost_fit = function(xgboost_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(xgboost_fit, newdata=xnew, type = "prob")
  list("prediction"=fit,"weights"=w)
}



model = xgboost_fit(as.matrix(X_training), as.matrix(W_training$W))
preds = predict.xgboost_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))

svm_fit = function(x,y,args=list()){
  model = do.call(svm, c(list(y = y, x = x, probability = TRUE, kernel = "sigmoid", 
                              type = 'C-classification', gamma = 0.1), args))
  return(model)
}

predict.svm_fit <- function(svm_fit, x, y, xnew = NULL, weights = FALSE) {
  if (is.null(xnew)) xnew = x
  if (weights) {
    w = NULL
  } else {
    w = NULL
  }
  
  # Perform predictions using the trained SVM model
  fit <- predict(svm_fit, newdata = xnew, probability = TRUE) %>% attr("probabilities")
  
  # Return the predictions and weights (if any)
  result <- list("prediction" = fit, "weights" = w)
  return(result)
}


model = svm_fit(x=as.matrix(X_training), y=as.matrix(W_training$W))
preds = predict.svm_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))





lda_fit = function(x,y,args=list()){
  model = do.call(lda, c(list(x = x, grouping = y, cv = TRUE)))
  model
}

predict.lda_fit = function(lda_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  fit <- predict(lda_fit, X_test, type="prob")$posterior
  list("prediction"=fit,"weights"=w)
  
}


qda_fit = function(x,y,args=list()){
  y = as.factor(y)
  
  model = do.call(qda, c(list(x = x, grouping = y, cv = TRUE), args))
  model
}

predict.qda_fit = function(qda_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit <- predict(qda_fit, X_test, type="prob")$posterior
  list("prediction"=fit,"weights"=w)
  
}


model <- qda_fit(x = as.matrix(X_training), y = as.factor(W_training$W))
preds <- predict.qda_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))


probability_forest_fit = function(x,y, args=list()){
  model = do.call(probability_forest, c(list(X = x, Y = y), args))
  model
}

predict.probability_forest_fit = function(probability_forest_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(probability_forest_fit, newdata=xnew)
  list("prediction"=fit,"weights"=w)
  
}



model <- probability_forest_fit(x = as.matrix(X_training), y = as.factor(W_training$W))
preds <- predict.probability_forest_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))


bagging_fit = function(x, y, args = list()) {
  data = as.data.frame(cbind(y, x))
  data$y <- as.factor(data$y)
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


predict.bagging_fit = function(bagging_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  data = as.data.frame(xnew)
  data$y <- as.factor(0)
  
  fit = predict(bagging_fit, newdata = data)$prob
  list("prediction"=fit,"weights"=w)
}


model = bagging_fit(x = as.matrix(X_training), y = W_training$W)
preds = predict.bagging_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))



boosting_fit = function(x, y, args = list()) {
  data = as.data.frame(cbind(y, x))
  data$y <- as.factor(data$y)
  model = do.call(boosting, c(list(formula = y ~ ., data = data,
                                  mfinal = 50,
                                  control = rpart.control(objective = "binary:logistic",
                                                          eval_metric = "rmse",
                                                          eta = 0.3, max_depth = 2,
                                                          minsplit = 10,  # Adjust minsplit to an appropriate value
                                                          cp = -1
                                  )), args))
  model
}

predict.boosting_fit = function(boosting_fit,x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  data = as.data.frame(xnew)
  data$y <- as.factor(0)
  
  fit = predict(boosting_fit, newdata = data)$prob
  list("prediction"=fit,"weights"=w)
}


model = boosting_fit(x = as.matrix(X_training), y = W_training$W)
preds = predict.boosting_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))


adaboost_fit = function(x,y,args=list()){
  data = as.data.frame(cbind(y, x))
  data$y <- as.factor(data$y)
  model = do.call(adaboost ,c(list(data=data, formula = y ~ ., nIter=20, loss="logistic", w=NULL), args))
  model
}

predict.adaboost_fit = function(adaboost_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(adaboost_fit, newdata=xnew, type = "probs")
  list("prediction"=fit,"weights"=w)
}


model = adaboost_fit(x = as.matrix(X_training), y = as.factor(W_training$W))
preds = predict.boosting_fit(model, x = as.matrix(X_training), y = as.matrix(W_training$W), xnew = as.matrix(X_test))












model_fit = function(x,y,args=list()){
  model = do.call( ... ,c(list(...), args))
  model
}

predict.model_fit = function(model_fit, x,y,xnew=NULL,weights=FALSE){
  if (is.null(xnew)) xnew = x
  if (weights==TRUE) {
    w = NULL
  }
  else w = NULL
  
  fit = predict(model_fit, )
  list("prediction"=fit,"weights"=w)
}