compute_mean_diff_RCT <- function(DF){
  RCT_ATE <- mean(DF[DF$A == 1 & DF$S == 1, "Y"]) - mean(DF[DF$A == 0  & DF$S == 1, "Y"])  
  return(RCT_ATE)
}

compute_ipsw <- function(DF, normalized = FALSE, estimation = "logit"){
  
  N <- nrow(DF)
  n <- nrow(DF[DF$S ==1, ])
  m <- nrow(DF[DF$S ==0, ])
  
  temp <- DF
  
  # Estimation of P(V = 1 | X)
  # p <-- P(V = 1 | X) 
  
  if (estimation == "logit"){
    
    # with logistic regression
    p.fit  <- glm(S ~., family = binomial("logit"), data = temp[, !names(temp) %in% c("A", "Y")])
    p <- predict(p.fit, type = "response", newdata = temp)
    
  } 
  else if (estimation == "forest") {
    break
  } 
  else {
    print("Estimation parameter should be forest or logit.")
    break
  }
  
  # Store odds
  temp$odds <- ((1 - p)/p)
  
  # Keep only RCT for the rest of the calculus
  temp <- temp[temp$S == 1,]
  
  if (normalized == FALSE){
    tau_ipsw <- (2/m)*with(temp, sum(odds*A*Y - odds*(1-A)*Y))  
  } else {
    tau_ipsw <- with(temp, sum(odds*A*Y/sum(odds*A) - odds*(1-A)*Y/sum(odds*(1-A))))
  }
  
  return(tau_ipsw)
}

# G-formula
compute_gformula <- function(DF, continuous_Y = TRUE){
  
  temp <- DF
  
  if(continuous_Y){
    mu_1 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 1, !names(temp) %in% c("S", "A")])
    mu_0 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 0, !names(temp) %in% c("S", "A")])
    
    mu_1_predict <- predict.lm(mu_1, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    mu_0_predict <- predict.lm(mu_0, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    
    tau_hat_gformula <- mean(mu_1_predict) - mean(mu_0_predict)
    
  } else 
    # binary outcome
    {
    mu_1 <- glm(Y ~., family = binomial("logit"), data = temp[temp$S == 1 & temp$A == 1, !names(temp) %in% c("S", "A")])
    mu_0 <- glm(Y ~., family = binomial("logit"), data = temp[temp$S == 1 & temp$A == 0, !names(temp) %in% c("S", "A")])
    
    mu_1_predict <- predict(mu_1, type = "response", newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    mu_0_predict <- predict(mu_0, type = "response", newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    
    tau_hat_gformula <- mean(mu_1_predict) - mean(mu_0_predict)
  }
  

  
  return(tau_hat_gformula)  
}


# AIPSW
compute_aipsw <- function(DF, normalized = FALSE) {
  
  N <- nrow(DF)
  n <- nrow(DF[DF$S ==1, ])
  m <- nrow(DF[DF$S ==0, ])
  
  temp <- DF
  
  ### Outcome modeling
  mu_1 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 1, !names(temp) %in% c("S", "A")])
  mu_0 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 0, !names(temp) %in% c("S", "A")])
  
  # RWE estimation
  mu_1_predict <- predict.lm(mu_1, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
  mu_0_predict <- predict.lm(mu_0, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
  tau_gformula <- mean(mu_1_predict) - mean(mu_0_predict)
  
  ### IPSW part
  # Estimation of P(V = 1 | X)
  # p <-- P(V = 1 | X) 
  # with logistic regression
  p.fit  <- glm(S ~., family = binomial("logit"), data = temp[, !names(temp) %in% c("A", "Y")])
  p <- predict(p.fit, type = "response", newdata = temp)
  
  # Store odds
  temp$odds <- ((1 - p)/p)
  
  #keep only rct
  temp <- temp[temp$V == 1,]
  
  temp$mu_11 <- predict.lm(mu_1, newdata = temp[, !names(temp) %in% c("S", "A")])
  temp$mu_10 <- predict.lm(mu_0, newdata = temp[, !names(temp) %in% c("S", "A")])
  
  if (normalized == FALSE){
    tau_ipsw <- (2/m)*with(temp, sum(odds*A*(Y - mu_11) - odds*(1-A)*(Y - mu_10)))  
  } else {
    tau_ipsw <- with(temp, sum(odds*A*(Y - mu_11)/sum(odds*A) - odds*(1-A)*(Y - mu_10)/sum(odds*(1-A))))
  }
  
  tau_aipsw <- tau_ipsw + tau_gformula
  return(tau_aipsw)
}

