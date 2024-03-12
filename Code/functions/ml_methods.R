
# multinomial logit glmnet
multinomial_logit <- list(type = "Classification",
                          library = "glmnet",
                          multiclass = TRUE) 
# define parameters
multinomial_logit$parameters <- data.frame(parameter = c("alpha", "lambda"),
                                           class = rep("numeric", 2),
                                           label = c("Mixing Percentage", "Regularization"))

multinomial_logit$fit <- logit_fit
multinomial_logit$predict <- predict.logit_fit
# design the parameter tuning grid
multinomial_logit$grid <- expand.grid(alpha = seq(0,1,0.05))#,


# multinomial logit nnet
multinomial_logit_nnet <- list(type = "Classification",
                               library = "nnet",
                               multiclass = TRUE) 
# define parameters
multinomial_logit_nnet$parameters <- data.frame(parameter = c("decay"),
                                                class = c("numeric"),
                                                label = c("weight decay"))

multinomial_logit_nnet$fit <- logit_nnet_fit
multinomial_logit_nnet$predict <- predict.logit_nnet_fit
# design the parameter tuning grid
multinomial_logit_nnet$grid <- expand.grid(decay = seq(0,1,0.05))


# KNN
knn <- list(type = "Classification",
            library = "kknn",
            multiclass = TRUE)
params <- data.frame(parameter = c("kmax", "distance", "kernel"),
                     class = c("numeric", "numeric", ""),
                     label = c("number of neighbors considered", 
                               "parameter of Minkowski distance",
                               "Kernel to use"))
knn$params <- params
knn$fit <- knn_fit
knn$predict <- predict.knn_fit
knn$grid <- expand.grid(kmax = seq(5,30,5),
                        #distance = seq(1,5,1),
                        kernel = c("rectangular", "gaussian", "optimal"))



# Naive Bayes bernulli
nb_bernulli <- list(
  type = "Classification",
  library = "naivebayes",
  multiclass = TRUE,
  params = data.frame(parameter = c("laplace"), 
                      class = c("numeric"), 
                      label = c("value of Laplace smoothing")),
  fit = nb_bernulli_fit,
  predict = predict.nb_bernulli_fit,
  grid = expand.grid(laplace = c(0, 1, 10, 50, 100, 200))
)

# Probability forest grf
probability_forest_grf <- list(
  type = "Classification",
  library = "grf",
  multiclass = TRUE,
  params = data.frame(parameter = c("mtry", "min.node.size", "alpha", "imbalance.penalty", "num.trees"), 
                      class = rep("numeric",5), 
                      label = c("Number of variables tried for each split",
                                "minimum number of observations in each tree leaf",
                                "controls maximum imbalance of split",
                                "controls how harshly imbalanced splits are penalized",
                                "number of trees grown"
                      )),
  fit = probability_forest_fit,
  predict = predict.probability_forest_fit,
  grid = expand.grid(min.node.size = seq(5, 20, 5),
                     alpha = seq(0,0.1,0.03),
                     imbalance.penalty = seq(0,1,0.5),
                     num.trees=1000
  )
)

# Probability forest ranger
probability_forest_ranger <- list(
  type = "Classification",
  library = "ranger",
  multiclass = TRUE,
  params = data.frame(parameter = c(""), class = c(""), label = c("")),
  fit = ranger_fit,
  predict = predict.ranger_fit,
  grid = expand.grid(
    mtry = ncol(X), #floor(x * c(.05, .15, .25, .333, .4)),
    min.node.size = seq(10, 20,5), 
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(.5, .66, .8, 1),
    splitrule = c("gini", "extratrees"),
    importance = c('none', 'impurity', 'permutation'),
    num.trees = 1000
  )
)

probability_forest_ranger$grid <- subset(probability_forest_ranger$grid, !(splitrule == "extratrees" & sample.fraction != 1))

# MLPC
mlpc <- list(
  type = "Classification",
  library = "nnet",
  multiclass = TRUE,
  params = data.frame(parameter = c("size"), class = c("numeric"), label = c("number of units in the hidden layer")),
  fit = mlpc_fit,
  predict = predict.mlpc_fit,
  grid = expand.grid(size=seq(1,10,1))
)

# xgboost
xgboost_method <- list(
  type = "Classification",
  library = "xgboost",
  multiclass = TRUE,
  params = data.frame(parameter = c("eta", "gamma", "nrounds", "min_child_weight", "max_depth", "subsample"), 
                      class = rep("numeric", 6), 
                      label = c("learning rate", "minimun loss reduction required for split", "# of boosting iterations", "minimum sum of instance weight (hessian) needed in a child", "maximum depth of a tree", "subsample ratio of the training instance")),
  fit = xgboost_fit,
  predict = predict.xgboost_fit,
  grid = expand.grid(eta = c(0.001,0.01,0.1,0.2),
                     gamma = seq(0,1,0.2),
                     nrounds = seq(5, 100, 5),
                     min_child_weight = seq(0,2,0.5),
                     max_depth = seq(3,10,2),
                     subsample = c(0.5,0.7,1))
)
# Randomly select 1000 rows
random_sample_indices <- sample(nrow(xgboost_method$grid), 200, replace = FALSE)
xgboost_method$grid <- xgboost_method$grid[random_sample_indices, ]
