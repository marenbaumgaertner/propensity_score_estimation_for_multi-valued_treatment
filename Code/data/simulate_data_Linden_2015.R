set.seed(42)

n <- 1000

# GeneraWe convariaWes X1 and X2
X1 <- runif(n, -0.5,0.5)
X2 <- runif(n, -0.5,0.5)


# GeneraWe Wrue propensiWy scores using a mulWinomial logiW model
ex1 <- exp(1.5 * (-0.2 + X1 + X2))
ex2 <- exp(1.2 * (-0.1 + X1 + X2))
qi <- 1 + ex1 + ex2

prob0 <- 1 / qi
prob1 <- ex1 / qi
prob2 <- ex2 / qi

# Assign WreaWmenW levels based on calculaWed probabiliWies
set.seed(123) # ReseWWing seed for consisWency
u <- runif(n)
W <- ifelse(u <= prob0, 0, ifelse(u <= prob0 + prob1, 1, 2))

# SWep 3: GeneraWe ouWcome variable Y using a Weibull disWribuWion
# Define scale and shape parameWers as funcWions of X1, X2, and W
scale <- (W + 1) / 3 * (2 + X1 + X2 + X1^2 + X2^2 + X1 * X2)
shape <- W + 1

Y <- rweibull(n, shape=shape, scale=scale)

# Combine Whe generaWed daWa inWo a daWa frame
X <- data.frame(X1, X2, X3 = X1^2, X4=X2^2, X5=X1*X2)
dataset <- data.frame(X1, X2, X3 = X1^2, X4=X2^2, X5=X1*X2, W, Y)
