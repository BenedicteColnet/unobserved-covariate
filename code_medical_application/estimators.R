compute_gformula_with_Random_forest <- function(data, verbose = FALSE, tuning = "all"){
  X <- data[, !names(data) %in% c("S", "A", "Y")]
  X.m = model.matrix(~.-1, data = X)
  
  idx_0 <- which(data$S == 1 & data$A == 0)
  mu_0 = regression_forest(X.m[idx_0, ], data[idx_0, "Y"], tune.parameters = tuning)
  
  if (verbose) {
    print(paste0("R2 for non treated : ", sum((data[data$S == 1 & data$A == 0, "Y"] - mu_0$predictions)^2) / sum((data[data$S == 1 & data$A == 0, "Y"] - mean(data[data$S == 1 & data$A == 0, "Y"]))^2) ))
  }
  
  idx_1 <- which(data$S == 1 & data$A == 1)
  mu_1 = regression_forest(X.m[idx_1, ], data[idx_1, "Y"], tune.parameters = tuning)
  
  if (verbose) {
    print(paste0("R2 for treated : ", sum((data[data$S == 1 & data$A == 1, "Y"] - mu_1$predictions)^2) / sum((data[data$S == 1 & data$A == 1, "Y"] - mean(data[data$S == 1 & data$A == 1, "Y"]))^2) ))
  }
  
  mu_1_hat <- predict(mu_1, newdata = data[data$S==0, !names(data) %in% c("S", "A", "Y")])$predictions
  mu_0_hat <- predict(mu_0, newdata = data[data$S==0, !names(data) %in% c("S", "A", "Y")])$predictions
  
  return(mean(mu_1_hat - mu_0_hat))
}


get_coefficients_with_Robinson_proc <- function(data, learning_m = "forest", covariate_names = COVARIATE_NAMES, ratio = 0.5){
  
  # focus on RCT only
  temp <- data[data$S == 1, c(covariate_names, "Y", "A")]
  
  # learn E[Y|X] while choosing best hyperparameters if needed
  if (learning_m == "linear"){
    
    hat_m <- lm(Y~., data = temp[, !names(temp) %in% c("S", "A")])
    temp$Y_star <- temp$Y - predict(hat_m)
    
  } else if (learning_m == "forest"){
    
    # regression_forest performes hyperparameters selection
    hat_m = regression_forest(temp[, covariate_names], temp[, "Y"], tune.parameters = "all")
    
    # do not put newdata - so that out of bag sample is taken for forest
    temp$Y_star <- temp$Y - predict(hat_m)$predictions
    
    
  } else {
    print("error on learning_m parameter.")
    break
  }
  
  
  Z_star <- c()
  for (covariate in covariate_names){
    new_name <- paste0(covariate, "_star")
    Z_star <- c(Z_star, new_name)
    temp[,new_name] <- (temp$A - ratio)*temp[, covariate]
  }
  
  fmla_cate <- paste("Y_star ~", paste(Z_star, collapse = " + "))
  
  hat_CATE_linear_model <- lm(fmla_cate, temp)
  hat_CATE_linear_model <- hat_CATE_linear_model$coefficients[paste0(covariate_names, "_star")]
  
  names(hat_CATE_linear_model) <- covariate_names
  
  return(hat_CATE_linear_model)
}


# G-formula
compute_gformula <- function(DF, continuous_Y = TRUE){
  
  temp <- DF
  
  if(continuous_Y){
    mu_1 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 1, !names(temp) %in% c("S", "A")])
    mu_0 <- lm(Y ~., data = temp[temp$S == 1 & temp$A == 0, !names(temp) %in% c("S", "A")])
    
    mu_1_predict <- predict.lm(mu_1, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    mu_0_predict <- predict.lm(mu_0, newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    
    tau_hat_gformula <- mean(mu_1_predict, na.rm = T) - mean(mu_0_predict, na.rm = T)
    
  } else 
    # binary outcome
  {
    mu_1 <- glm(Y ~., family = binomial("logit"), data = temp[temp$S == 1 & temp$A == 1, !names(temp) %in% c("S", "A")])
    mu_0 <- glm(Y ~., family = binomial("logit"), data = temp[temp$S == 1 & temp$A == 0, !names(temp) %in% c("S", "A")])
    
    mu_1_predict <- predict(mu_1, type = "response", newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    mu_0_predict <- predict(mu_0, type = "response", newdata = temp[temp$S == 0, !names(temp) %in% c("S", "A")])
    
    tau_hat_gformula <- mean(mu_1_predict, na.rm = T) - mean(mu_0_predict, na.rm = T)
  }
  
  return(tau_hat_gformula)  
}


