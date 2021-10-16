## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)
# Set random generator seed for reproducible results
set.seed(123)
# Libraries
library(ggplot2) # plots
library(tidyverse)
library(lubridate) # date management
library(caret)
library(rpart)
library(grf)

# Load methods from estimators.R 
source("./estimators.R")


## -----------------------------------------------------------------------------
star <- read.csv("./STAR_data.csv")
star$g1tchid <- as.factor(star$g1tchid)


## -----------------------------------------------------------------------------
star$Y = star$g1tlistss + star$g1treadss + star$g1tmathss
star$Y <- star$Y/3
star <- star %>% drop_na("Y")
star <- star %>% drop_na("g1classtype")
star <- star[star$g1classtype != 3,]
star$g1classtype <- ifelse(star$g1classtype == 2, 0, 1)
Kallus_X <- c("gender", "race", "birthmonth", "birthday", "birthyear", "gkfreelunch", "g1tchid", "g1freelunch", "g1surban")
star <- star %>% mutate_at(Kallus_X, ~replace_na(., 0))


## -----------------------------------------------------------------------------
assertthat::are_equal(nrow(star[star$g1classtype == 0,]), 2413)
assertthat::are_equal(nrow(star[star$g1classtype == 1,]), 1805)
assertthat::are_equal(nrow(star), 4218)
assertthat::are_equal(mean(star$Y), 540.0987829935199)
N = nrow(star)
ggplot(star, aes(x = Y, color = as.factor(g1classtype))) +
  geom_histogram(bins = 50, fill = "white") +
  facet_grid(~g1classtype) +
  theme_bw()


## -----------------------------------------------------------------------------
star$date <- paste(star$birthmonth, star$birthday, star$birthyear, sep = "-")
star$date  <- mdy(star$date)
star$age <- interval(star$date, "1985-01-01") %/% months(1)


## -----------------------------------------------------------------------------
names(star)[names(star) == "g1classtype"] <- "A"


## -----------------------------------------------------------------------------
X <- c("gender", "race", "age", "gkfreelunch", "g1tchid", "g1freelunch", "g1surban")
X_without_g1surban <- c("gender", "race", "age", "gkfreelunch", "g1tchid", "g1freelunch")
covariates <- c(X, "A", "Y")
star <- star[, covariates]
star$g1surban <- as.numeric(star$g1surban)


## -----------------------------------------------------------------------------
# Filter treatment / control observations, pulls outcome variable as a vector
y_treated <- star[star$A == 1, "Y"] 
y_control <- star[star$A == 0, "Y"] 
  
n_treated <- nrow(star[star$A == 1,]) 
n_control <- nrow(star[star$A == 0,]) 
  
# Difference in means
estimated_difference <- mean(y_treated) - mean(y_control)
  
# 95% Confidence intervals
se_hat <- sqrt(var(y_treated)/(n_treated-1) + var(y_control)/(n_control-1) )
lower_ci <- estimated_difference - 1.96 * se_hat
upper_ci <- estimated_difference + 1.96 * se_hat
  
GROUND_TRUTH = c(Difference = estimated_difference, CI_inf = lower_ci, CI_sup = upper_ci)
GROUND_TRUTH


## -----------------------------------------------------------------------------
RWE_index <- sample(1:N, 500, replace = FALSE)
not_RWE_index <- setdiff(1:N, RWE_index)
RWE <- star[RWE_index,]
star <- star[not_RWE_index,]

assertthat::are_equal(nrow(star)+nrow(RWE), N)

# control that the sample is a rather good estimate
t.test(RWE[RWE$A == 1, "Y"], RWE[RWE$A == 0, "Y"])

RWE$S <- rep(0, nrow(RWE))


## -----------------------------------------------------------------------------
class <- star %>% 
  group_by(g1tchid) %>%
  summarise(mean_g1surban = mean(g1surban), mean_race = mean(race), mean_lunch = mean(g1freelunch))
etas <- as.vector(as.matrix(class[, c("mean_g1surban", "mean_race", "mean_lunch")]) %*% c(-1, 0, 0)) + 0.5
ps = 1 / (1 + exp(-etas)) # P(S = 1 | X)
RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
class$S <- RCT_indicator
class$S <- ifelse(is.na(class$S), 0, class$S)
star <- merge(star, class[, c("S", "g1tchid")], by = "g1tchid")
RCT <- star[star$S == 1,]


## -----------------------------------------------------------------------------
ggplot(RCT, aes(x = race)) +
  geom_bar()
ggplot(RWE, aes(x = race)) +
  geom_bar()
ggplot(RCT, aes(x = g1freelunch)) +
  geom_bar()
ggplot(RWE, aes(x = g1freelunch)) +
  geom_bar()
ggplot(RCT, aes(x = g1surban)) +
  geom_bar()
ggplot(RWE, aes(x = g1surban)) +
  geom_bar()


## -----------------------------------------------------------------------------
star <- rbind(RCT, RWE)
star <- na.omit(star)


## -----------------------------------------------------------------------------
# Filter treatment / control observations, pulls outcome variable as a vector
y_treated_RCT <- RCT[RCT$A == 1, "Y"] 
y_control_RCT <- RCT[RCT$A == 0, "Y"] 
  
n_treated_RCT <- nrow(RCT[RCT$A == 1,]) 
n_control_RCT <- nrow(RCT[RCT$A == 0,]) 
  
# Difference in means
estimated_difference <- mean(y_treated_RCT) - mean(y_control_RCT)
  
# 95% Confidence intervals
se_hat <- sqrt(var(y_treated_RCT)/(n_treated_RCT-1) + var(y_control_RCT)/(n_control_RCT-1) )
lower_ci <- estimated_difference - 1.96 * se_hat
upper_ci <- estimated_difference + 1.96 * se_hat
  
TAU1 = c(Difference = estimated_difference, CI_inf = lower_ci, CI_sup = upper_ci)
TAU1


## ---- echo=T------------------------------------------------------------------
stratified_bootstrap <- function(DF, nboot = 500){
  
  estimands <- c()
  
  for (i in 1:nboot){
  
    # random resamples from RCT
    n = nrow(DF[DF$S == 1,])
    index_RCT = sample(1:n, n, replace = TRUE)
    
    # random resamples from RWD
    m = nrow(DF[DF$S == 0,])
    index_RWD = sample(1:m, m, replace = TRUE)
    
    # new data set
    RCT_RWD <- rbind(DF[which(DF$S==1),][index_RCT,],
                     DF[which(DF$S==0),][index_RWD,])
    # estimation
    estimands <- c(estimands, compute_gformula_with_Random_forest(RCT_RWD))
  }
  return(estimands)
}

ate <- compute_gformula_with_Random_forest(star[, c("gender", "race", "g1freelunch", "gkfreelunch", "g1surban", "S", "Y", "A")])
bootstrap <- stratified_bootstrap(DF=star[, c("gender", "race", "g1freelunch", "gkfreelunch", "g1surban", "S", "Y", "A")])
bootstrap_sorted <- sort(bootstrap)
CI_inf <- quantile(bootstrap_sorted, 0.025)
CI_sup <- quantile(bootstrap_sorted, 0.975)
GENERALIZED_ATE_all_cov = c(Difference = ate, CI_inf = CI_inf, CI_sup = CI_sup)
GENERALIZED_ATE_all_cov

ate <- compute_gformula_with_Random_forest(star[, c("gender", "race", "g1freelunch", "gkfreelunch", "S", "Y", "A")])
bootstrap <- stratified_bootstrap(DF=star[, c("gender", "race", "g1freelunch", "gkfreelunch",  "S", "Y", "A")])
bootstrap_sorted <- sort(bootstrap)
CI_inf <- quantile(bootstrap_sorted, 0.025)
CI_sup <- quantile(bootstrap_sorted, 0.975)
GENERALIZED_ATE_miss_g1surban = c(Difference = ate, CI_inf = CI_inf, CI_sup = CI_sup)
GENERALIZED_ATE_miss_g1surban


ate <- compute_gformula_with_Random_forest(star[, c("gender", "race", "S", "Y", "A")])
bootstrap <- stratified_bootstrap(DF=star[, c("gender", "race", "S", "Y", "A")])
bootstrap_sorted <- sort(bootstrap)
CI_inf <- quantile(bootstrap_sorted, 0.025)
CI_sup <- quantile(bootstrap_sorted, 0.975)
GENERALIZED_ATE_only_gender = c(Difference = ate, CI_inf = CI_inf, CI_sup = CI_sup)
GENERALIZED_ATE_only_gender


## -----------------------------------------------------------------------------
results <- rbind(GROUND_TRUTH, TAU1, GENERALIZED_ATE_all_cov, GENERALIZED_ATE_miss_g1surban)
results <- data.frame(names = row.names(results), results)
results$names <- factor(results$names,
                        levels = c("GROUND_TRUTH", "TAU1", "GENERALIZED_ATE_all_cov", "GENERALIZED_ATE_miss_g1surban"))


## -----------------------------------------------------------------------------
results$names_recoded = case_when(results$names == "GROUND_TRUTH" ~ "STAR RCT",
                                  results$names == "TAU1" ~ "Biaised RCT",
                                  results$names == "GENERALIZED_ATE_all_cov" ~"Generalized ATE \n (all covariates)",
                                  results$names == "GENERALIZED_ATE_miss_g1surban" ~ "Generalized ATE \n (without g1surban)")

ggplot(results, aes(x=names_recoded, y=Difference)) + 
  geom_point(size = 2)+
  geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=.2,
                 position=position_dodge(0.05)) +
  geom_hline(yintercept = GROUND_TRUTH["Difference"], linetype = "dashed", color = "red") + 
  coord_flip() +
  theme_bw() +
  xlab("") +
  ylab("Estimated ATE") +
  theme(text = element_text(size=13))
  ggsave("./figures/bilan_star_with_random_forest.pdf", width = 4, height = 2.5)

