P = 5
SIGMA = matrix(c(1, 0, 0, 0, 0.8,
                 0, 1, 0, 0, 0,
                 0, 0, 1, 0, 0,
                 0, 0, 0, 1, 0,
                 0.8, 0, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
DELTA = c(30, 30, -10, 0, 0)
BETA_s = c(-0.4, 0, -0.3, -0.3, 0) 
BETA = c(5, 5, 5, 5, 5)
COVARIATE_NAMES = c("X1", "X2", "X3", "X4", "X5")
TRUE_ATE = sum(DELTA)