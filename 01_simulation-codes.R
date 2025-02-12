###############################################################################
## Project: Collapsibility 
## Program 1: create and run simulation function 

# Goal:
#### create simulation function 
#### save output from simulations 

## Created by: Ya-Hui Yu (Jan 29, 2025)
# This program is modified from Ashley's dgm exploration program 
################################################################################

#===============================================================================
# iteration=1
# simulate data from linear/logistic models
simulation_function <- function(iteration, add_interaction = TRUE){
  
#-------------------------------------------------------------------------------
#### DGM 
  n <- 1e3
  W <- runif(n)
  A <- rbinom(n, 1, plogis(W))
  
  # Generate outcomes based on whether interaction term is added or not
  if (add_interaction) {
    # Include interaction term A*W
    logit_Y <- rbinom(n, 1, plogis(A + W + A * W))
    continuous_Y <- A + W + A * W + rnorm(n)
  } else {
    # Without interaction term
    logit_Y <- rbinom(n, 1, plogis(A + W))
    continuous_Y <- A + W + rnorm(n)
  }
  
  linear_Y <- as.numeric(continuous_Y > median(continuous_Y))
  data <- data.frame(monte_carlo = iteration, W, A, logit_Y, linear_Y)
  
#-------------------------------------------------------------------------------  
#### g-computation 
  
  analytic_function_gcomp <- function(outcome_model = model1){
    # these are used throughout to compute g-comp
    data_Ais1 <- data; data_Ais1$A <- 1
    data_Ais0 <- data; data_Ais0$A <- 0  
    
    pred_Ais1 <- predict(
      outcome_model, newdata = data_Ais1, type = "response"
    )
    pred_Ais0 <- predict(
      outcome_model, newdata = data_Ais0, type = "response"
    )
    
    # E[Y(1)] estimate (for ATE)
    EY1_hat <- mean(pred_Ais1)
    # E[Y(0)] estimate (for ATT)
    EY0_hat <- mean(pred_Ais0)
    # E[Y(1) | A = 1] estimate (for ATT)
    EY1_Ais1_hat <- mean(pred_Ais1[A == 1])
    # E[Y(0) | A = 1] estimate (for ATT)
    EY0_Ais1_hat <- mean(pred_Ais0[A == 1])
    
    # contrasts 1
    ATE_hat <- EY1_hat - EY0_hat
    ATT_hat <- EY1_Ais1_hat - EY0_Ais1_hat
    
    # contrasts 2
    ATE_logOR_hat <- log((EY1_hat/(1 - EY1_hat))/(EY0_hat/(1 - EY0_hat)))
    ATT_logOR_hat <- log((EY1_Ais1_hat/(1 - EY1_Ais1_hat))/
                           (EY0_Ais1_hat/(1 - EY0_Ais1_hat)))
    
    res <- data.frame(ATE_hat = ATE_hat,
                      ATT_hat = ATT_hat,
                      ATE_logOR_hat = ATE_logOR_hat,
                      ATT_logOR_hat = ATT_logOR_hat,
                      comparisonRD = ATE_hat - ATT_hat,
                      comparisonOR = ATE_logOR_hat - ATT_logOR_hat)
    return(res)
  }
  
  #.............................................................................
  #### Different scenarios for fitting models
  
  # case 1, linear DGM, single linear regression, no interactions
  model1 <- glm(
    linear_Y ~ A + W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 2, linear DGM, single linear regression, interactions
  model2 <- glm(
    linear_Y ~ A*W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 3, logit DGM, single linear regression, no interactions
  model3 <- glm(
    logit_Y ~ A+W, 
    data = data,
    family = gaussian(link = "identity")
  )
  # case 4, logit DGM, single linear regression, interactions
  model4 <- glm(
    logit_Y ~ A*W, 
    data = data,
    family = gaussian(link = "identity")
  )
  
  # case 5, linear DGM, single logit regression, no interactions
  model5 <- glm(
    linear_Y ~ A + W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 6, linear DGM, single logit regression, interactions
  model6 <- glm(
    linear_Y ~ A*W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 7, logit DGM, single logit regression, no interactions
  model7 <- glm(
    logit_Y ~ A+W, 
    data = data,
    family = binomial(link = "logit")
  )
  # case 8, logit DGM, single logit regression, interactions
  model8 <- glm(
    logit_Y ~ A*W, 
    data = data,
    family = binomial(link = "logit")
  )
  
  model_list_gcomp <- list(model1, model2, model3, model4, model5, 
                           model6, model7, model8)
  #.............................................................................
  #### obtain final results from gcomp (sim_res_gcomp)
  simulation_results <- 
    lapply(model_list_gcomp, 
           function (x) analytic_function_gcomp(outcome_model = x))
  
  sim_res_gcomp <- do.call(rbind, simulation_results)
  sim_res_gcomp$scenario <- 1:8
  sim_res_gcomp$method <- "gcomp"
  
#------------------------------------------------------------------------------- 
#### ipw
  analytic_function_ipw <- function(outcome, linkfunction){
  
  data$pscore <- glm(A ~ W, data = data, 
                     family = binomial(link = "logit"))$fitted.values
  
  #marginal probability of A 
  pa<- mean(data$A)
  data$sw <- (pa/data$pscore)*data$A +
    (mean(1 - data$A)/(1 - data$pscore))*(1 - data$A)
  
  data$sw_att <- ifelse(data$A == 1, 1, 
                        pa / (1 - pa) * (1 - data$pscore) / data$pscore)
  
  ATE_hat <- summary(
    glm(outcome ~ A,
        data = data,
        weights = sw,
        family = binomial(link = linkfunction )))$coefficients[2,1]

  ATT_hat <- summary(
   glm(outcome ~ A,
      data = data,
      weights = sw_att,
      family = binomial(link = linkfunction)))$coefficients[2,1]
  
  estimate <- ATE_hat - ATT_hat
  
  results <- c(ATE_hat, ATT_hat, estimate) 
  return(results)
  
  }
 
  #.............................................................................  
  #### Different scenarios for fitting models
  ## QUESTIONS: cannot have both RD and OR for each scenario; can we include A*W?
  ## (these estimates were set to 9999 as n/a)
  
  # case 1, linear DGM, single linear regression, no interactions
  r1<-analytic_function_ipw(outcome=linear_Y, linkfunction="identity")
  s1<-data.frame(ATE_hat = r1[[1]],
                 ATT_hat = r1[[2]],
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = r1[[3]],
                 comparisonOR = 9999,
                 scenario=1,
                 method="ipw")

  # case 2, linear DGM, single linear regression, interactions (X)
  s2<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = 9999,
                 comparisonOR = 9999,        
                 scenario=2,
                 method="ipw")

  # case 3, logit DGM, single linear regression, no interactions
  r3<-analytic_function_ipw(outcome=logit_Y, linkfunction="identity")
  s3<-data.frame(ATE_hat = r3[[1]],
                 ATT_hat = r3[[2]],
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = r1[[3]],
                 comparisonOR = 9999,
                 scenario=3,
                 method="ipw")
  
  # case 4, logit DGM, single linear regression, interactions (X)
  s4<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = 9999,
                 comparisonOR = 9999,        
                 scenario=4,
                 method="ipw")
  
  # case 5, linear DGM, single logit regression, no interactions
  r5<-analytic_function_ipw(outcome=linear_Y, linkfunction="logit")
  s5<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = r5[[1]],
                 ATT_logOR_hat = r5[[2]],
                 comparisonRD = 9999,
                 comparisonOR = r5[[3]],
                 scenario=5,
                 method="ipw")
  
  # case 6, linear DGM, single logit regression, interactions (X)
  s6<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = 9999,
                 comparisonOR = 9999,        
                 scenario=6,
                 method="ipw")
  
  # case 7, logit DGM, single logit regression, no interactions
  r7<-analytic_function_ipw(outcome=logit_Y, linkfunction="logit")
  s7<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = r7[[1]],
                 ATT_logOR_hat = r7[[2]],
                 comparisonRD = 9999,
                 comparisonOR = r7[[3]],
                 scenario=7,
                 method="ipw")
  
  # case 8, logit DGM, single logit regression, interactions (X)
  s8<-data.frame(ATE_hat = 9999,
                 ATT_hat = 9999,
                 ATE_logOR_hat = 9999,
                 ATT_logOR_hat = 9999,
                 comparisonRD = 9999,
                 comparisonOR = 9999,        
                 scenario=8,
                 method="ipw")
  #.............................................................................
  #### obtain final results from gcomp (sim_res_ipw)
  sim_res_ipw <- rbind(s1, s2, s3, s4, s5, s6, s7, s8)

#------------------------------------------------------------------------------- 
### AIPW and TMLE (these two methods use same set of predicted values 
### from same set of SL)
  
  analytic_function_aipwtmle <- function(outcome, sl.lib){

  Y <- outcome
  A <- data$A
  covariates_ps <- data %>% select (W)
  covariates_mu <- data %>% select (A,W)
  
  bound_func <- function(probs, lower_bound, upper_bound){
    probs <- if_else(probs < lower_bound, lower_bound, probs)
    probs <- if_else(probs > upper_bound, upper_bound, probs)
    return(probs)
  }
  
  num.folds <- 10
  n<-nrow(data)
  folds <- sort(seq(n) %% num.folds) + 1
  fold_dat <- tibble(id = 1:n, folds)
  fold_index <- split(fold_dat$id, fold_dat$folds)
  
  #create a wrap for "SL.lm.interaction" (todo)
  
  fit_mu <- CV.SuperLearner(Y = Y,
                            X = covariates_mu, 
                            method = "method.NNLS", 
                            family = binomial(),
                            SL.library = sl.lib,
                            cvControl = list(V = num.folds, 
                                             validRows = fold_index),
                            control = list(saveCVFitLibrary = F),
                            parallel = "seq",
                            verbose = T)
  
  fit_pi <- CV.SuperLearner(Y = A,
                            X = covariates_ps,
                            method = "method.NNLS", 
                            family = binomial(),
                            SL.library = sl.lib,
                            cvControl = list(V = num.folds, 
                                             validRows = fold_index),
                            control = list(saveCVFitLibrary = F),
                            parallel = "seq",
                            verbose = T)
  
  pscore <- as.matrix(fit_pi$SL.predict)
  pscore <- bound_func(pscore, .025, .975)

  # output coefficients from SL outcome model
  coef_SL_mu <- NULL
  for(i in 1:num.folds){
    coef_SL_mu <- rbind(coef_SL_mu, fit_mu$AllSL[[i]]$coef)
  }
    
  mu_hat1 <- NULL
  for(i in 1:num.folds){
    covs_temp <- covariates_mu[fold_index[[1]], ] %>%mutate(!!"A" := 1)
    mu_hat1 <- rbind(mu_hat1, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = covs_temp, 
                             onlySL=T)$pred)
  }
  
  mu_hat0 <- NULL
  for(i in 1:num.folds){
    covs_temp <- covariates_mu[fold_index[[1]], ] %>%mutate(!!"A" := 0)
    mu_hat0 <- rbind(mu_hat0, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = covs_temp, 
                             onlySL=T)$pred)
  }
  
  mu_hat1 <- bound_func(mu_hat1, .0125, .975)
  mu_hat0 <- bound_func(mu_hat0, .0125, .975)
  
  #............................................................................. 
  #### AIPW 
  mu1_aipw <- (A/ pscore) * (outcome - mu_hat1) + mu_hat1
  mu0_aipw <- ((1 - A) / (1 - pscore)) * (outcome - mu_hat0) + mu_hat0
  
  W<- data$W
  aipw<- data.frame(mu1_aipw, mu0_aipw, pscore, mu_hat1, mu_hat0, W, A, Y)

  # contrasts 1
  ATE_hat <- mean(mu1_aipw - mu0_aipw)
  ATT_hat <- mean(aipw[aipw$A == 0, ]$mu1_aipw) - 
                mean (aipw[aipw$A == 0, ]$mu0_aipw)
  
  # contrasts 2
  ATE_logOR_hat <-log((mean(mu1_aipw) / mean(mu0_aipw)) / 
                        ((1 - mean(mu1_aipw)) / (1 - mean(mu0_aipw))))
  ATT_logOR_hat <- 
    log((mean(aipw[aipw$A == 1, ]$mu1_aipw) / mean(aipw[aipw$A == 1, ]$mu0_aipw)) / 
    ((1 - mean(aipw[aipw$A == 1, ]$mu1_aipw)) / (1 - mean(aipw[aipw$A == 1, ]$mu0_aipw))))
 
  
  res_aipw <- data.frame(ATE_hat = ATE_hat,
                         ATT_hat = ATT_hat,
                         ATE_logOR_hat = ATE_logOR_hat,
                         ATT_logOR_hat = ATT_logOR_hat,
                         comparisonRD = ATE_hat - ATT_hat,
                         comparisonOR = ATE_logOR_hat - ATT_logOR_hat,
                         method = "AIPW") 
  
  #.............................................................................
  #### TMLE (we will use the predicted valued from SL when we estimate aipw,
  # therefore, we will use aipw dataset here)
  
  tmle_est <- tmle(Y = aipw$Y, 
                   A = aipw$A, 
                   W = as.matrix(aipw$W),
                   family = "binomial",
                   Q =cbind(aipw$mu_hat0, aipw$mu_hat1),
                   g1W = pscore,
                   verbose = T)
  # contrasts 1
  ATE_hat <- tmle_est$estimates$ATE$psi
  ATT_hat <- tmle_est$estimates$ATT$psi
  
  # contrasts 2
  ATE_logOR_hat <- tmle_est$estimates$OR$log.psi
   
  res_tmle <- data.frame(ATE_hat = ATE_hat,
                         ATT_hat = ATT_hat,
                         ATE_logOR_hat = ATE_logOR_hat,
                         ATT_logOR_hat = 9999,
                         comparisonRD = ATE_hat - ATT_hat,
                         comparisonOR = 9999,
                         method = "TMLE") 

  res_both <-rbind(res_aipw, res_tmle)
  results <-list(res_both, coef_SL_mu)

  return(results)
  
  }

  #.............................................................................  
  #### Different scenarios for fitting models
  
  # case 1, linear DGM, single linear regression, no interactions
  s1<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.lm"))[[1]]
  
  # case 2, linear DGM, single linear regression, interactions 
  # (X- need sl.lm.interaction)
  s2<-data.frame(ATE_hat = c(9999,9999),
                 ATT_hat = c(9999,9999),
                 ATE_logOR_hat = c(9999,9999),
                 ATT_logOR_hat = c(9999,9999),
                 comparisonRD = c(9999,9999),
                 comparisonOR = c(9999,9999),        
                 method= c("aipw","tmle"))
  
  # case 3, logit DGM, single linear regression, no interactions
  s3<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.lm"))[[1]]
  
  # case 4, logit DGM, single linear regression, interactions 
  # (X- need sl.lm.interaction)
  s4<-data.frame(ATE_hat = c(9999,9999),
                 ATT_hat = c(9999,9999),
                 ATE_logOR_hat = c(9999,9999),
                 ATT_logOR_hat = c(9999,9999),
                 comparisonRD = c(9999,9999),
                 comparisonOR = c(9999,9999),        
                 method= c("aipw","tmle"))
  
  # case 5, linear DGM, single logit regression, no interactions
  s5<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.glm"))[[1]]
  
  # case 6, linear DGM, single logit regression, interactions 
  s6<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.glm","SL.glm.interaction"))[[1]]
  
  # case 7, logit DGM, single logit regression, no interactions
  s7<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.glm"))[[1]]
  
  # case 8, logit DGM, single logit regression, interactions 
  s8<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.glm","SL.glm.interaction"))[[1]]
  
  # case 9, linear DGM, sl.lib with both algorithms (need to save coef of sl)
  s9<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.lm",
                                          "SL.glm"))[[1]]
  
  s9_c<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                   sl.lib=c("SL.lm",
                                            "SL.glm"))[[2]]
  
  # case 10, logit DGM, sl.lib with both algorithms (need to save coef of sl)
  s10<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                  sl.lib=c("SL.lm",
                                           "SL.glm"))[[1]]
  
  s10_c<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                    sl.lib=c("SL.lm",
                                             "SL.glm"))[[2]]
  
  # case 11, linear DGM, sl.lib with all algorithms (need to save coef of sl)
  s11<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.lm",
                                          "SL.glm","SL.glm.interaction"))[[1]]
  
  s11_c<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.lm",
                                          "SL.glm","SL.glm.interaction"))[[2]]
  
  # case 12, logit DGM, sl.lib with all algorithms (need to save coef of sl)
  s12<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.lm",
                                          "SL.glm","SL.glm.interaction"))[[1]]
  
  s12_c<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                  sl.lib=c("SL.lm",
                                           "SL.glm","SL.glm.interaction"))[[2]]

  #.............................................................................
  #### obtain final results from AIPW and TMLE (sim_res_dr)
  sim_res_dr<-rbind(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12)
  sim_res_dr$ scenario <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
                            11,11,12,12)
  
  coef_SL2_mu<- data.frame(rbind(s9_c, s10_c))
  coef_SL2_mu$scenario <- rep(c(9, 10), each = 10)
  coef_SL3_mu<- data.frame(rbind(s11_c, s12_c))
  coef_SL3_mu$scenario <- rep(c(11, 12), each = 10)
#-------------------------------------------------------------------------------
#### obtain final results for this function
  
  table<-cbind(iteration, rbind(sim_res_gcomp, sim_res_ipw, sim_res_dr))
  coef2_table <-cbind(iteration, coef_SL2_mu)
  coef3_table <-cbind(iteration, coef_SL3_mu)
  final_results_list<-list (table, coef2_table,coef3_table)
  
  return(final_results_list)
}



#===============================================================================
### Run simulation for DGM with interaction
set.seed(123)

# sim_res <- mclapply(1:4, function(x) simulation_function(iteration = x,
#                                                          add_interaction = FALSE), 
#                     mc.cores = 3, 
#                     mc.preschedule = FALSE, 
#                     mc.cleanup = TRUE)

sim_res <- lapply(1:200, function(x) simulation_function(iteration = x, 
                                                       add_interaction = TRUE))
res_inter <- do.call(rbind, lapply(sim_res, `[[`, 1))
coef2_inter<- do.call(rbind, lapply(sim_res, `[[`, 2))
coef3_inter<- do.call(rbind, lapply(sim_res, `[[`, 3))

write_csv(coef2_inter, here("output", "coef2_inter.csv"))
write_csv(coef3_inter, here("output", "coef3_inter.csv"))

#...............................................................................
# Results for RD,OR comparison
results <- res_inter %>% 
  group_by(scenario,method) %>% 
  summarize(RD_diff = mean(comparisonRD),
            OR_diff = mean(comparisonOR),
            ATE_rd = mean(ATE_hat),
            ATE_or = mean(ATE_logOR_hat),
            ATT_rd = mean(ATT_hat),
            ATT_or = mean(ATT_logOR_hat))

#results[ ]<- lapply(results, function(x) if(is.numeric(x)) round(x, 4) else x)

# Define the scenarios 
scenario_table <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  AnalyticModel = c("linear", "Linear", "Linear", "Linear", 
                    "Logit", "Logit", "Logit", "Logit",
                    "SL2","SL2","SL3","SL3"),
  InteractionInModel = c("No", "Yes", "No", "Yes", 
                         "No", "Yes", "No", "Yes",
                         "No","No","Yes","Yes")
)

results_inter <- scenario_table %>%
            right_join(results, by = "scenario") %>%
            mutate_if(is.numeric, ~ na_if(., 9999))

write_csv(results_inter, here("output", "results_inter.csv"))

#-------------------------------------------------------------------------------
### Run simulation for DGM without interaction
set.seed(123)
# sim_res <- mclapply(1:1, function(x) simulation_function(iteration = x,
#                                                          add_interaction = FALSE), 
#                     mc.cores = 3, 
#                     mc.preschedule = FALSE, 
#                     mc.cleanup = TRUE)

sim_res_no <- lapply(1:200, function(x) simulation_function(iteration = x, 
                                                       add_interaction = FALSE))
res_nointer <- do.call(rbind, lapply(sim_res_no, `[[`, 1))
coef2_nointer<- do.call(rbind, lapply(sim_res_no, `[[`, 2))
coef3_nointer<- do.call(rbind, lapply(sim_res_no, `[[`, 3))

write_csv(coef2_nointer, here("output", "coef2_nointer.csv"))
write_csv(coef3_nointer, here("output", "coef3_nointer.csv"))

#...............................................................................
# Results for RD,OR comparison
results <- res_nointer %>% 
  group_by(scenario,method) %>% 
  summarize(RD_diff = mean(comparisonRD),
            OR_diff = mean(comparisonOR),
            ATE_rd = mean(ATE_hat),
            ATE_or = mean(ATE_logOR_hat),
            ATT_rd = mean(ATT_hat),
            ATT_or = mean(ATT_logOR_hat))

#results[ ]<- lapply(results, function(x) if(is.numeric(x)) round(x, 4) else x)

# Define the scenarios [Confirm with each model]
scenario_table <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  AnalyticModel = c("linear", "Linear", "Linear", "Linear", 
                    "Logit", "Logit", "Logit", "Logit",
                    "SL2","SL2","SL3","SL3"),
  InteractionInModel = c("No", "Yes", "No", "Yes", 
                         "No", "Yes", "No", "Yes",
                         "No","No","Yes","Yes")
)

results_nointer <- scenario_table %>%
  right_join(results, by = "scenario") %>%
  mutate_if(is.numeric, ~ na_if(., 9999))

write_csv(results_nointer, here("output", "results_nointer.csv"))


#-------------------------------------------------------------------------------
## Combine two tables 
results_inter <- cbind(DGM_interaction = "Yes", results_inter)
results_nointer <- cbind(DGM_interaction = "No", results_nointer)
final_table <- rbind(results_inter,results_nointer)

write_csv(final_table, here("output", "results_final.csv"))

