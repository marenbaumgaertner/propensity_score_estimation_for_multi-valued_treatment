#set.seed(42)

n <- 2000

# Generate convariates X1 and X2
X1 <- runif(n, -0.5,0.5)
X2 <- runif(n, -0.5,0.5)


# Generate true propensity scores using a multinomial logit model
ex1 <- exp(1.5 * (-0.2 + X1 + X2))
ex2 <- exp(1.2 * (-0.1 + X1 + X2))
qi <- 1 + ex1 + ex2

prob0 <- 1 / qi
prob1 <- ex1 / qi
prob2 <- ex2 / qi

# Assign treatment levels based on calculated probabilities
u <- runif(n, 0, 1)
W <- ifelse(u <= prob0, 0, ifelse(u <= prob0 + prob1, 1, 2))

# Generate outcome variable Y using a Weibull distribution
# Define scale and shape parameters as functions of X1, X2, and W
scale <- (W + 1) / 3 * (2 + X1 + X2 + X1^2 + X2^2 + X1 * X2)
shape <- W + 1

Y <- rweibull(n, shape=shape, scale=scale)

# Combine the generated data into a data frame
X <- data.frame(X1, X2, X3 = X1^2, X4=X2^2, X5=X1*X2)
dataset <- data.frame(X1, X2, X3 = X1^2, X4=X2^2, X5=X1*X2, W, Y)
X=as.matrix(X)
table(dataset$W)
