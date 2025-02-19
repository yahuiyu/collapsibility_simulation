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
iteration=1
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
  
  #create a wrap for "SL.lm.interaction
  SL.lm.interaction <- 
    function (Y, X, newX, family, obsWeights, model = TRUE, ...)  {
    if (is.matrix(X)) {
      X = as.data.frame(X)
    }
    fit <- stats::lm(Y ~ .^2, data = X, weights = obsWeights, model = model)
    if (is.matrix(newX)) {
      newX = as.data.frame(newX)
    }
    pred <- predict(fit, newdata = newX, type = "response")
    if (family$family == "binomial") {
      pred = pmin(pmax(pred, 0), 1)
    }
    fit <- list(object = fit, family = family)
    class(fit) <- "SL.lm"
    out <- list(pred = pred, fit = fit)
    return(out)
  }
  
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
                            SL.library = "SL.glm",
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
  
  # mu_hat1 <- bound_func(mu_hat1, .0125, .975)
  # mu_hat0 <- bound_func(mu_hat0, .0125, .975)
  
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
  
  # Extract TMLE estimates to calculate ATT_logOR_hat
  Q1 <- tmle_est$Qstar[,2]  
  Q0 <- tmle_est$Qstar[,1]  
  dat_attor <- data.frame(cbind(A, Q1, Q0))
  
  odds1_ATT <-mean(dat_attor$Q1[dat_attor$A==1]) /
              (1-mean(dat_attor$Q1[dat_attor$A==1]))
  
  odds0_ATT <-mean(dat_attor$Q0[dat_attor$A==0]) /
              (1-mean(dat_attor$Q0[dat_attor$A==0]))
  
  ATT_logOR_hat <- log(odds1_ATT/odds0_ATT)
   
  res_tmle <- data.frame(ATE_hat = ATE_hat,
                         ATT_hat = ATT_hat,
                         ATE_logOR_hat = ATE_logOR_hat,
                         ATT_logOR_hat = ATT_logOR_hat,
                         comparisonRD = ATE_hat - ATT_hat,
                         comparisonOR = ATE_logOR_hat - ATT_logOR_hat,
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
  s2<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                      sl.lib=c("SL.lm","SL.lm.interaction"))[[1]]
  
  # case 3, logit DGM, single linear regression, no interactions
  s3<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.lm"))[[1]]
  
  # case 4, logit DGM, single linear regression, interactions 
  s4<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.lm","SL.lm.interaction"))[[1]]
  
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
                                 sl.lib=c("SL.lm", "SL.lm.interaction",
                                          "SL.glm","SL.glm.interaction"))[[1]]
  
  s11_c<-analytic_function_aipwtmle(outcome=data$linear_Y, 
                                 sl.lib=c("SL.lm", "SL.lm.interaction",
                                          "SL.glm","SL.glm.interaction"))[[2]]
  
  # case 12, logit DGM, sl.lib with all algorithms (need to save coef of sl)
  s12<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                 sl.lib=c("SL.lm", "SL.lm.interaction",
                                          "SL.glm","SL.glm.interaction"))[[1]]
  
  s12_c<-analytic_function_aipwtmle(outcome=data$logit_Y, 
                                  sl.lib=c("SL.lm", "SL.lm.interaction",
                                           "SL.glm","SL.glm.interaction"))[[2]]

  #.............................................................................
  #### obtain final results from AIPW and TMLE (sim_res_dr)
  sim_res_dr<-rbind(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12)
  sim_res_dr$ scenario <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
                            11,11,12,12)
  
  coef_SL2_mu<- data.frame(rbind(s9_c, s10_c))
  coef_SL2_mu$scenario <- rep(c(9, 10), each = 10)
  coef_SL4_mu<- data.frame(rbind(s11_c, s12_c))
  coef_SL4_mu$scenario <- rep(c(11, 12), each = 10)
#-------------------------------------------------------------------------------
#### obtain final results for this function
  
  table<-cbind(iteration, rbind(sim_res_gcomp, sim_res_dr))
  coef2_table <-cbind(iteration, coef_SL2_mu)
  coef4_table <-cbind(iteration, coef_SL4_mu)
  final_results_list<-list (table, coef2_table,coef4_table)
  
  return(final_results_list)
}



#===============================================================================
### Run simulation for DGM with interaction
set.seed(123)
sim_res <- mclapply(1:500, function(x) simulation_function(iteration = x,
                                                         add_interaction = TRUE), 
                   mc.cores = 3, 
                   mc.preschedule = FALSE, 
                   mc.cleanup = TRUE)

# sim_res <- lapply(1:2, function(x) simulation_function(iteration = x, 
#                                                        add_interaction = TRUE))

res_inter <- do.call(rbind, lapply(sim_res, `[[`, 1))
coef2_inter<- do.call(rbind, lapply(sim_res, `[[`, 2))
coef4_inter<- do.call(rbind, lapply(sim_res, `[[`, 3))

current_date <- Sys.Date()
filename <- here("output", paste0("res_inter_full", current_date, ".csv"))
write_csv(res_inter, filename)

current_date <- Sys.Date()
filename <- here("output", paste0("coef2_inter_", current_date, ".csv"))
write_csv(coef2_inter, filename)

current_date <- Sys.Date()
filename <- here("output", paste0("coef4_inter_", current_date, ".csv"))
write_csv(coef4_inter, filename)

write_csv(coef2_inter, here("output", "coef2_inter.csv"))

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

# Define the scenarios 
scenario_table <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  AnalyticModel = c("linear", "Linear", "Linear", "Linear", 
                    "Logit", "Logit", "Logit", "Logit",
                    "SL2","SL2","SL4","SL4"),
  InteractionInModel = c("No", "Yes", "No", "Yes", 
                         "No", "Yes", "No", "Yes",
                         "No","No","Yes","Yes")
)

results_inter <- scenario_table %>%
            right_join(results, by = "scenario") %>%
            mutate_if(is.numeric, ~ na_if(., 9999))

current_date <- Sys.Date()
filename <- here("output", paste0("results_inter_", current_date, ".csv"))
write_csv(results_inter, filename)

#-------------------------------------------------------------------------------
### Run simulation for DGM without interaction
set.seed(123)
sim_res_no <- mclapply(1:500, function(x) 
  simulation_function(iteration = x, add_interaction = FALSE), 
  mc.cores = 3, 
                   mc.preschedule = FALSE, 
                   mc.cleanup = TRUE)

# sim_res_no <- lapply(1:2, function(x) simulation_function(iteration = x, 
#                                                        add_interaction = FALSE))

res_nointer <- do.call(rbind, lapply(sim_res_no, `[[`, 1))
coef2_nointer<- do.call(rbind, lapply(sim_res_no, `[[`, 2))
coef4_nointer<- do.call(rbind, lapply(sim_res_no, `[[`, 3))

current_date <- Sys.Date()
filename <- here("output", paste0("res_nointer_full", current_date, ".csv"))
write_csv(res_nointer, filename)

current_date <- Sys.Date()
filename <- here("output", paste0("coef2_nointer_", current_date, ".csv"))
write_csv(coef2_nointer, filename)

current_date <- Sys.Date()
filename <- here("output", paste0("coef4_nointer_", current_date, ".csv"))
write_csv(coef4_nointer, filename)


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

# Define the scenarios [Confirm with each model]
scenario_table <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  AnalyticModel = c("linear", "Linear", "Linear", "Linear", 
                    "Logit", "Logit", "Logit", "Logit",
                    "SL2","SL2","SL4","SL4"),
  InteractionInModel = c("No", "Yes", "No", "Yes", 
                         "No", "Yes", "No", "Yes",
                         "No","No","Yes","Yes")
)

results_nointer <- scenario_table %>%
  right_join(results, by = "scenario") %>%
  mutate_if(is.numeric, ~ na_if(., 9999))


current_date <- Sys.Date()
filename <- here("output", paste0("results_nointer_", current_date, ".csv"))
write_csv(results_nointer, filename)

#-------------------------------------------------------------------------------
## Combine two tables 
results_inter <- cbind(DGM_interaction = "Yes", results_inter)
results_nointer <- cbind(DGM_interaction = "No", results_nointer)
final_table <- rbind(results_inter,results_nointer)

final_table <- final_table %>% 
  select(scenario, DGM, DGM_interaction, AnalyticModel, InteractionInModel,
         method, RD_diff, OR_diff, ATE_rd, ATT_rd, ATE_or, ATT_or) 

current_date <- Sys.Date()
filename <- here("output", paste0("results_final_", current_date, ".csv"))
write_csv(final_table, filename)

#===============================================================================
### True values 
set.seed(123)
n <- 5e7
W <- runif(n)
A <- rbinom(n, 1, plogis(W))

logit_Y <- rbinom(n, 1, plogis(A + W + A * W))
continuous_Y <- A + W + A * W + rnorm(n)

# With interaction
continous_Y_inter <- A + W +  A*W + rnorm(n)
threshold_inter <- median(continous_Y_inter)

EY1_logit_inter <- rbinom(n, 1, plogis(1 + W + 1 * W))
continuous_Y1_inter <- 1 + W + 1 * W + rnorm(n)
EY1_linear_inter <- as.numeric(continuous_Y1_inter > 
                                   median(threshold_inter))

EY0_logit_inter <-  logit_Y <- rbinom(n, 1, plogis(0 + W + 0 * W))
continuous_Y0_inter <- 0 + W + 0 * W + rnorm(n)
EY0_linear_inter <- as.numeric(continuous_Y0_inter > 
                                   median(threshold_inter))

# Without interaction
continous_Y_nointer <- A + W + rnorm(n)
threshold_nointer <- median(continous_Y_nointer)

EY1_logit_nointer <- rbinom(n, 1, plogis(1 + W))
continuous_Y1_nointer <- 1 + W + rnorm(n)
EY1_linear_nointer <- as.numeric(continuous_Y1_nointer >     
                                     median(threshold_nointer))
                                   
EY0_logit_nointer <- rbinom(n, 1, plogis(0 + W))
continuous_Y0_nointer <- 0 + W + rnorm(n)
EY0_linear_nointer <- as.numeric(continuous_Y0_nointer >     
                                     median(threshold_nointer))

#...............................................................................
ATE_linear_inter <-  mean(EY1_linear_inter)- mean(EY0_linear_inter)  
ATE_logit_inter <- mean(EY1_logit_inter) - mean(EY0_logit_inter)  

ATE_linear_nointer <- mean(EY1_linear_nointer) - mean(EY0_linear_nointer)  
ATE_logit_nointer <- mean(EY1_logit_nointer) - mean(EY0_logit_nointer)  
  
ATT_linear_inter <- mean(EY1_linear_inter[A == 1]) - mean(EY0_linear_inter[A == 1])  
ATT_logit_inter <- mean(EY1_logit_inter[A == 1]) - mean(EY0_logit_inter[A == 1])    

ATT_linear_nointer <- mean(EY1_linear_nointer[A == 1]) - mean(EY0_linear_nointer[A == 1])  
ATT_logit_nointer <- mean(EY1_logit_nointer[A == 1]) - mean(EY0_logit_nointer[A == 1])    

RD_diff_logit_inter <- ATE_logit_inter-ATT_logit_inter
RD_diff_logit_nointer <- ATE_logit_nointer-ATT_logit_nointer

RD_diff_linear_inter <- ATE_linear_inter-ATT_linear_inter
RD_diff_linear_nointer <- ATE_linear_nointer-ATT_linear_nointer

#-------------------------------------------------------------------------------
# Create the table based on the results 

# Calculate bias for ATE and ATT
res_nointer_full <- read_csv(here("output","res_nointer_full2025-02-15.csv"))
res_inter_full <- read_csv(here("output","res_inter_full2025-02-15.csv"))
res_final<- read_csv(here("output","results_final_2025-02-15.csv"))

# Define the scenarios [Confirm with each model]
scenario_table_nointer <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  True_ATE = 
    c(ATE_linear_nointer, ATE_linear_nointer, ATE_logit_nointer, ATE_logit_nointer,
      ATE_linear_nointer, ATE_linear_nointer, ATE_logit_nointer, ATE_logit_nointer,
      ATE_linear_nointer, ATE_logit_nointer, ATE_linear_nointer, ATE_logit_nointer),
  
  True_ATT = 
    c(ATT_linear_nointer, ATT_linear_nointer, ATT_logit_nointer, ATT_logit_nointer,
      ATT_linear_nointer, ATT_linear_nointer, ATT_logit_nointer, ATT_logit_nointer,
      ATT_linear_nointer, ATT_logit_nointer, ATT_linear_nointer, ATT_logit_nointer)
)

scenario_table_inter <- data.frame(
  scenario = 1:12,
  DGM = c("Linear", "Linear", "Logit", "Logit", 
          "Linear", "Linear", "Logit", "Logit",
          "Linear", "Logit","Linear", "Logit"),
  True_ATE = c(ATE_linear_inter, ATE_linear_inter, ATE_logit_inter, ATE_logit_inter,
               ATE_linear_inter, ATE_linear_inter, ATE_logit_inter, ATE_logit_inter,
               ATE_linear_inter, ATE_logit_inter, ATE_linear_inter, ATE_logit_inter),
  
  True_ATT = c(ATT_linear_inter, ATT_linear_inter, ATT_logit_inter, ATT_logit_inter,
               ATT_linear_inter, ATT_linear_inter, ATT_logit_inter, ATT_logit_inter,
               ATT_linear_inter, ATT_logit_inter, ATT_linear_inter, ATT_logit_inter)
)

results_full_nointer <- scenario_table_nointer %>%
  right_join(res_nointer_full, by = "scenario") 
results_full_nointer$ATE_bias <- 
  results_full_nointer$ATE_hat -results_full_nointer$True_ATE
results_full_nointer$ATT_bias <- 
  results_full_nointer$ATT_hat -results_full_nointer$True_ATT

results_full_inter <- scenario_table_inter %>%
  right_join(res_inter_full, by = "scenario") 
results_full_inter$ATE_bias <- 
  results_full_inter$ATE_hat -results_full_inter$True_ATE
results_full_inter$ATT_bias <- 
  results_full_inter$ATT_hat -results_full_inter$True_ATT

results_full_nointer1 <- results_full_nointer %>% 
  group_by(scenario,method) %>% 
  summarize(ATE_rdbias = mean(ATE_bias),
            ATT_rdbias = mean(ATT_bias))

results_full_inter1 <- results_full_inter %>% 
  group_by(scenario,method) %>% 
  summarize(ATE_rdbias = mean(ATE_bias),
            ATT_rdbias = mean(ATT_bias))

results_full_nointer1 <- cbind(DGM_interaction = "No", results_full_nointer1)
results_full_inter1 <- cbind(DGM_interaction = "Yes", results_full_inter1)
final_table1 <- rbind(results_full_inter1,results_full_nointer1)
results_final_withbias <- final_table1 %>%
  right_join(res_final, by = c("scenario","DGM_interaction","method"))

results_final_withbias <- results_final_withbias %>% 
  select(scenario, DGM, DGM_interaction, AnalyticModel, InteractionInModel,
         method, RD_diff, OR_diff, ATE_rdbias,ATT_rdbias,
         ATE_rd, ATT_rd, ATE_or, ATT_or) 

write_csv(results_final_withbias,
          here("output","result_final_withbias.csv"))

res_final<- read_csv(here("output","results_final_2025-02-15.csv"))



