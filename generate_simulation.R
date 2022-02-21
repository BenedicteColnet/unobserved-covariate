# note that correlation is let as a parameters because useful for some simulations
# by default the CATE is linear followint the parameters given in BETA and DELTA
generate_simulation <- function(bs = BETA_s, beta = BETA, delta = DELTA, m = 10000, correlation = 0.8, model = "linear", ratio = 0.5){
  
  # covariates definition
  p = length(bs)
  X_names <- paste("X", 1:p, sep = "")
  covariates_names <- c(X_names)
  Sigma = matrix(c(1, 0, 0, 0, correlation,
                   0, 1, 0, 0, 0,
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0,
                   correlation, 0, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
  
  mu = rep(1, p)
  
  # generate source population for RCT
  source_data_for_RCT <- mvrnorm(n = m, mu = mu, Sigma = Sigma)  
  source_data_for_RCT <- as.data.frame(source_data_for_RCT)
  names(source_data_for_RCT) <- covariates_names
  
  # sample RCT
  etas <- as.vector(as.matrix(source_data_for_RCT[, paste("X", 1:p, sep = "")]) %*% bs)
  
  # P(S = 1 | X) 
  ps = 1 / (1 + exp(-etas))
  source_data_for_RCT$ps <- ps
  RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
  source_data_for_RCT$S <- RCT_indicator
  
  # random treatment assignment within the RCT
  source_data_for_RCT$A <- ifelse(source_data_for_RCT$S == 1, rbinom(nrow(source_data_for_RCT), 1, ratio), NA)
  
  # keep only interesting variables
  source_data_for_RCT <- source_data_for_RCT[, c(covariates_names, "A", "S")]
  
  # drop other data
  RCT <- source_data_for_RCT[source_data_for_RCT$S == 1,]
  rm(source_data_for_RCT)
  
  # generate target population
  RWD <-  mvrnorm(n = m, mu = mu, Sigma = Sigma)
  RWD <- as.data.frame(RWD)
  names(RWD) <- covariates_names
  
  RWD$S <- rep(0, m)
  RWD$A <- rep(NA, m)
  
  # stack RCT and RWE
  DF <- rbind(RCT, RWD)
  
  # reset row number
  rownames(DF) <- 1:nrow(DF)
  
  # generate nonparametric function
  
  if (model == "linear"){
    DF$Y = beta[1]*DF$X1 + beta[2]*DF$X2 + beta[3]*DF$X3 + beta[4]*DF$X4 + beta[5]*DF$X5 
    
  } else if (model == "pmax"){
    DF$Y = pmax(0, DF$X1, DF$X2) + 2 * (DF$X3 - 0.5)^2 + DF$X4 + 0.5 * DF$X5
    
  } else if (model == "sin"){
    DF$Y = 10*sin(DF$X1*DF$X2)
    
  } else if (model == "log"){
    DF$Y = 10*log(abs(DF$X1)) +DF$X2*DF$X4 + 1/(exp(abs(DF$X3)))
    
  } else {
    stop("Error: model parameters does not exist.") }
  
  # add CATE and error
  error = rnorm(n = nrow(DF), mean = 0, sd = 0.2)
  DF$Y =  DF$Y +
    delta[1]*(DF$A == 1)*DF$X1 +
    delta[2]*(DF$A == 1)*DF$X2 +
    delta[3]*(DF$A == 1)*DF$X3 +
    delta[4]*(DF$A == 1)*DF$X4 +
    delta[5]*(DF$A == 1)*DF$X5 +
    error
  
  return(DF)  
}


# generate totally non-parametric simulation
# X1 and X2 are treatment effect modifier
# X1 and X4 are linked to trial inclusion
generate_totally_nonparametric_simulation <- function(m = 50000, p = 5, ratio = 0.5){
  
  # covariates definition
  X_names <- paste("X", 1:p, sep = "")
  covariates_names <- c(X_names)
  Sigma = matrix(c(1, 0, 0, 0, 0,
                   0, 1, 0, 0, 0,
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), nrow = 5, ncol = 5, byrow = TRUE)
  
  mu = rep(1,p)
  
  # generate source population for RCT
  source_data_for_RCT <- mvrnorm(n = m, mu = mu, Sigma = Sigma)  
  source_data_for_RCT <- as.data.frame(source_data_for_RCT)
  names(source_data_for_RCT) <- covariates_names
  
  # sample RCT
  ps <- rep(0.05, nrow(source_data_for_RCT))
  source_data_for_RCT$ps <- ps
  source_data_for_RCT$ps <- ifelse(source_data_for_RCT$X1 > 1 & source_data_for_RCT$X4 > 1, 0.4,
                                   ifelse(source_data_for_RCT$X1 > 1 & source_data_for_RCT$X4 <= 1, 0.5, 
                                          ifelse(source_data_for_RCT$X3 > 1, 0.35, 0.05)))
  
  ps <- source_data_for_RCT$ps
  RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
  source_data_for_RCT$S <- RCT_indicator
  
  # random treatment assignment within the RCT
  source_data_for_RCT$A <- ifelse(source_data_for_RCT$S == 1, rbinom(nrow(source_data_for_RCT), 1, ratio), NA)
  
  # keep only interesting variables
  source_data_for_RCT <- source_data_for_RCT[, c(covariates_names, "A", "S")]
  
  # drop other data
  RCT <- source_data_for_RCT[source_data_for_RCT$S == 1,]
  rm(source_data_for_RCT)
  
  # generate target population
  RWD <-  mvrnorm(n = m, mu = mu, Sigma = Sigma)
  RWD <- as.data.frame(RWD)
  names(RWD) <- covariates_names
  
  RWD$S <- rep(0, m)
  RWD$A <- rep(NA, m)
  
  # stack RCT and RWE
  DF <- rbind(RCT, RWD)
  
  # reset row number
  rownames(DF) <- 1:nrow(DF)
  
  # add CATE and error (X1 and X2 are treatment effect modifier)
  error = rnorm(n = nrow(DF), mean = 0, sd = 0.1)
  
  DF$Y_0 <- 5*DF$X1*DF$X1 +  4*DF$X3 + 3*pmax(1,DF$X5) + log(1 + 3*abs(DF$X2)) + error
  
  DF$Y_1 <- DF$Y_0 +  3*exp(-DF$X4) + DF$X3*DF$X3 + 10*DF$X1
  
  DF$Y <- ifelse(DF$A == 1, DF$Y_1, DF$Y_0)
  
  return(DF)  
}

