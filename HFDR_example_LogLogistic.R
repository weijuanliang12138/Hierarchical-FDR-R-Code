#
# Hierarchical False Discovery Rate Control for High-dimensional Survival Analysis with Interactions
# Authors: Weijuan Liang, Qingzhao Zhang, Shuangge Ma
# Date: Oct 8th, 2023
#

#-------------------------- Source packages and functions ----------------------
# Source the debiased estimator packages
library(flare)
library(ncvreg)
library(MASS)
library(plyr)
library(magrittr)
source("lasso_inference.R")
# Source other required packages

# task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# seed    <- task_id
# print(seed)
# set.seed(seed)


#-------------------------- initialization parameters --------------------------
# initialization all parameters 
alpha_set   <- c(0.01, 0.05, 0.10, 0.15, 0.2, 0.25, 0.30)
censor_rate <- 0.2
rho         <- 0.3
n           <- 500
S           <- 10
A           <- 1
d           <- 200
q           <- 5


#--------------------------------- Functions -----------------------------------
# Construct the Kaplan-Meier weights
Weight_func <- function(delta_sort, n){
  w         = vector()
  w[1]      = delta_sort[1]/n
  for(i in 2:n) w[i] = delta_sort[i]/(n-i+1)*prod(sapply(c(1:(i-1)), function(j){((n-j)/(n-j+1))^{delta_sort[j]}}))
  return(w)
}

# Critical value for the non-hierary FCD-Surv
critical_t_original <- function(t, p, z){
  2*p*(1-pnorm(t))/(max(sum(ifelse(abs(z) >= t, 1, 0)), 1))
}

# Critical value for the proposed method
critical_t_proposed <- function(t, z, d = d, q = q){
  R_t_main_alpha = which(abs(z[1:d]) >= t)
  R_t_main_gamma = which(abs(z[(d+1):(d+q)]) >= t)
  
  R_t_interaction_indx = vector()
  for(i in R_t_main_alpha){
    for(j in 1:q){
      R_t_interaction_indx = c(R_t_interaction_indx, abs(z[d+i*q+j]) > t)
    }
  }
  R_t = length(R_t_main_alpha) + length(R_t_main_gamma) + sum(R_t_interaction_indx)
  (d*2*(1-pnorm(t))*(1+q*2*(1-pnorm(t))))/max(R_t, 1)
}

P_value_func <- function(Phi_star_j, Y_star){
  fit_j      = lm(Y_star ~ Phi_star_j)
  summary(fit_j)$coefficients[2,4]
}

#--------------------------------- Generate data -------------------------------
# Generate the data
## Generate the X, with AR correlation
Mean  <- rep(0, d)
Sigma <- rho^{abs(outer(1:d, 1:d, "-"))}
X     <- mvrnorm(n, Mean, Sigma)

## Generate Z
Z     <- replicate(q, rnorm(n))

## Construct the Phi
Phi                <- matrix(nrow = n, ncol = p)
Phi[, 1:d]         <- X
Phi[, (d+1):(d+q)] <- Z
for(i in 1:n) Phi[i,(d+q+1):p] = kronecker(X[i,], Z[i,])

## Define the true parameters
alpha_true_index <- 1:S
gamma_true_index <- c(2,5)
beta_true_index  <- rep(alpha_true_index*q, each = 2)+d+rep(gamma_true_index, times = S)
index_true       <- c(alpha_true_index, gamma_true_index + d, beta_true_index)
true_theta       <- rep(0, p)
true_theta[index_true[1:S]]    <- 2
true_theta[d+gamma_true_index] <- 2
true_theta[beta_true_index]    <- 1

# Generate the event time and censoring time
# Log-Logistic
u_0   <- runif(n = n)
TT    <- u_0/(1-u_0)*exp(A*(Phi %*% true_theta))
CC    <- quantile(TT, 1-censor_rate)
## observed Y
Y     <- ifelse(log(TT) <= log(CC), log(TT), log(CC))
## censor indicator
delta <- ifelse(TT <= CC, 1, 0)

#------------------------------- KM Weighted design-----------------------------
KM_weight_func <- function(Y, Phi, delta){
  ## Sort the data according to Y
  observe      = data.frame(Y, Phi, delta)
  observe_sort = arrange(observe, Y)
  ## subtract the sorted data
  Y_sort       = observe_sort$Y
  Phi_sort     = as.matrix(observe_sort[,2:(p+1)], nrow = n, ncol = p)
  delta_sort   = observe_sort$delta
  ## Calculate the weight function
  W            = Weight_func(delta_sort, n)
  W_half       = diag(sqrt(n*W))
  ## The weighted Design matrix and response
  Phi_star     = crossprod(W_half, as.matrix(Phi_sort))
  Y_star       = crossprod(W_half, Y_sort)
  
  ret_list     = list(Y_sort = Y_sort,Y_star = Y_star, Phi_sort =  Phi_sort, 
                      Phi_star = Phi_star, delta_sort = delta_sort)
  return(ret_list)
}


#-------------------------- Test Statistics by Debiased Estimator --------------
Debiased_Est_Test_Stat_func  <- function(Y_sort, Y_star, Phi_sort, Phi_star, delta_sort){
  #------------------ Fit debiased Lasso ------------------
  fit          = SSLasso(Phi_star, Y_star, mu = 2*sqrt(log(p)/n), verbose = T)
  theta        = fit$coef
  
  #---------------- Calculate covariance matrix -----------
  ## hat(phi)
  phi_hat      = Phi_sort * as.vector(Y_sort - Phi_sort %*% theta)
  
  ## hat(tau)_0
  tau_hat_0    = vector(length = n)
  for(i in 1:n){
    index_k      = which((Y_sort < Y_sort[i])&(delta_sort == 0))
    sum1         = 0
    for(k in index_k){
      sum1.      = sum1 + 1/(n-sum(Y_sort <= Y_sort[k]))}
    tau_hat_0[i] = exp(sum1)
  }
  
  ## hat(tau)_1
  tau_hat_1   = matrix(nrow = n, ncol = p)
  for(i in 1:n){
    index_k2      = which((Y_sort > Y_sort[i])&(delta_sort == 1))
    sum1          = 0
    for(k in index_k2){
      sum1        = sum1 + phi_hat[k,] * tau_hat_0[k]/(n-sum(Y_sort <= Y_sort[i]))} 
    tau_hat_1[i,] = sum1
  }
  
  ## hat(tau)_2
  tau_hat_2   = matrix(nrow = n, ncol = p)
  for(i in 1:n){
    index_k3   = which((Y_sort < Y_sort[i])&(delta_sort == 0))
    sum1       = 0
    for(k in index_k3){
      index_l  = which((Y_sort > Y_sort[k])&(delta_sort == 1))
      numert   =  apply(phi_hat[index_l,]*tau_hat_0[index_l], 2, sum)
      denomnt  = (n-sum(Y_sort <= Y_sort[k]))^2
      sum1     = sum1 + numert/denomnt
    }
    tau_hat_2[i,] = sum1
  }
  
  ## zeta_hat
  zeta_hat   = phi_hat * tau_hat_0 * delta_sort + tau_hat_1 * (1-delta_sort) - tau_hat_2
  
  ## covariance matrix
  Sigma_hat  = cov(zeta_hat, zeta_hat)
  M_hat      = fit$M[-1,-1]
  Lambad_hat = M_hat %*% Sigma_hat %*% t(M_hat)
  
  #----------------- Construct test statistics ------------
  theta_thrshed = theta
  theta_thrshed[which(theta <= 1e-3)] = 0
  
  Test_stat  = sqrt(n) * theta_thrshed / sqrt(diag(Lambad_hat))
  ret_list   = list(theta_thrshed = theta_thrshed, Test_stat = Test_stat)
  return(ret_list)
}


BH_p_value_func <- function(Y_star, Phi_star, alpha, index_true, hierarchy = "FALSE"){
  #------------------ Marginal p-value -----------------
  if(hierarchy == "FALSE"){
    P_val_adj    = p.adjust(apply(Phi_star, 2, P_value_func, Y_star), method = "BH")
    S_m1      = which(P_val_adj < alpha)
    TP_m1     = length(which(S_m1 %in% index_true))
    FP_m1     = length(S_m1) - TP_m1
    Miss_m1   = length(index_true) - length(which(index_true %in% S_m1))  
    ## FDR
    FDR_m1    = FP_m1/length(S_m1)
    Power_m1  = TP_m1/length(index_true)
    
    ret_list  = list(P_val_adj = P_val_adj, S_m1 = S_m1, TP_m1 = TP_m1,
                     FP_m1 = FP_m1, Miss_m1 = Miss_m1, FDR_m1 = FDR_m1,
                     Power_m1 = Power_m1)
  }else if(hierarchy == "TRUE"){
    # Level 1 test for alpha and gamma 
    P_val_main  = p.adjust(apply(Phi_star[, 1:(d+q)], 2, P_value_func, Y_star), method = "BH")
    S_alpha     = which(P_val_main[1:d] < alpha)
    S_gamma     = which(P_val_main[1:q+d] < alpha) + d
    S_interac   = rep(d+S_alpha*q, each = q) + 1:q
    
    if(length(S_interac) != 0){
      P_val_interac = p.adjust(apply(Phi_star[, S_interac], 2, P_value_func, Y_star), method = "BH")
    }else{P_val_interac = 1}
    
    S_beta      = S_interac[P_val_interac < alpha]
    
    ## Results
    S_m2     = c(S_alpha, S_gamma, S_beta)
    TP_m2    = length(which(S_m2 %in% index_true))
    FP_m2    = length(S_m2) - TP_m2
    Miss_m2  = length(index_true) - length(which(index_true %in% S_m2))
    ## FDR
    FDR_m2   = FP_m2/length(S_m2)
    Power_m2 = TP_m2/length(index_true)
    
    ret_list  = list(P_val_main = P_val_main, P_val_interac = P_val_interac,
                     S_m2 = S_m2, TP_m2 = TP_m2,
                     FP_m2 = FP_m2, Miss_m2 = Miss_m2, FDR_m2 = FDR_m2,
                     Power_m2 = Power_m2)
  }
  return(ret_list)
}


#-----------------------------------  FDR control ------------------------------
#' @param alpha, prespecified FDR level
#' @param n, sample size
#' @param S, number of nonzero main genes
#' @param rho, correlation among covariats X
#' @param A, scalar 
#' @param d, number of total main genes
#' @param q, number of total environment covariates
#' @param censor_rate, censoring rate of survival
#' @param Y, survival outcome
#' @param Phi, augment covariates
#' @param delta, censoring indicator
#' 
#' @return FDR, empirical FDR
#' @return Power, empirical Power
#' @return Miss, number of missed true important covariates
#' @return TP, true positive
#' @return FP, false positive
#' @return MSE, mean square errors of all covariates

LogLogistic_fun <- function(alpha, n, S, rho, A, d, q, censor_rate, Y, Phi, delta){
  # Manipulate the data
    data_list   = KM_weight_func(Y, Phi, delta)
    Y_sort      = data_list$Y_sort
    Y_star      = data_list$Y_star
    Phi_sort    = data_list$Phi_sort
    Phi_star    = data_list$Phi_star
    delta_sort  = data_list$delta_sort
  
  #------------------- marginal methods -----------------------
  # Method 1 BH (without hierarchy)
    ret1      = BH_p_value_func(Y_star = Y_star, Phi_star = Phi_star, alpha = alpha,
                                index_true = index_true, hierarchy = "FALSE")
    S_m1      = ret1$S_m1
    TP_m1     = ret1$TP_m1
    FP_m1     = ret1$FP_m1
    Miss_m1   = ret1$Miss_m1
    ## FDR
    FDR_m1    = ret1$FDR_m1
    Power_m1  = ret1$Power_m1
  
  # Method 2 BH-hierarchy (hierarchy)
    ret2      = BH_p_value_func(Y_star = Y_star, Phi_star = Phi_star, alpha = alpha, 
                                index_true = index_true, hierarchy = "TRUE")
    S_m2      = ret2$S_m2
    TP_m2     = ret2$TP_m2
    FP_m2     = ret2$FP_m2
    Miss_m2   = ret2$Miss_m2
    ## FDR
    FDR_m2    = ret2$FDR_m2
    Power_m2  = ret2$Power_m2
  
  #------------------- Debiased methods -----------------------  
    Est_Test_Stat = Debiased_Est_Test_Stat_func(Y_sort = Y_sort, Y_star = Y_star, 
                                                Phi_sort = Phi_sort, Phi_star = Phi_star, 
                                                delta_sort = delta_sort)
    # Test statistics and estimated theta
    Test_stat     = Est_Test_Stat$Test_stat
    theta_thrshed = Est_Test_Stat$theta_thrshed
    
    tp            = (2*log(p) - 2*log(log(p)))^{1/2}
    t_step        = seq(0, tp, length = 1000)
  
  # Method 3 Proposed (hierarchy)
    # Level 1 test for alpha and beta
    t01       = try(t_step[min(which(sapply(t_step, critical_t_proposed, z = Test_stat, d = d, q = q) < alpha))])
    if("try-error" %in% class(t01)){t01 = sqrt(2*log(p))}
    
    # Select results
    S_alpha   = which(abs(Test_stat[1:d]) > t01)
    S_gamma   = which(abs(Test_stat[(d+1):(d+q)]) > t01) + d
    S_interac = rep(d+S_alpha*q, each = q) + 1:q
    S_beta    = S_interac[abs(Test_stat[S_interac]) > t01]
    
    # Results
    S_m3      = c(S_alpha, S_gamma, S_beta)
    TP_m3     = length(which(S_m3 %in% index_true))
    FP_m3     = length(S_m3) - TP_m3
    Miss_m3   = length(index_true) - length(which(index_true %in% S_m3))
    FDR_m3    = FP_m3/length(S_m3)
    Power_m3  = TP_m3/length(index_true)
    
    # MSE
    theta_propsd        = theta_thrshed
    theta_propsd[-S_m3] = 0
    MSE_propsd          = mean((theta_propsd - true_theta)^2)
  
  # Method 4 Debiased method (without hierarchy)
    t02      = try(t_step[min(which(sapply(t_step, critical_t_original, p = p, z = Test_stat) < alpha))])
    if("try-error" %in% class(t02) || is.na(t02)){t02 = sqrt(2*log(p))}
    
    S_m4     = which(Test_stat > t02)
    TP_m4    = length(which(S_m4 %in% index_true)) 
    FP_m4    = length(S_m4) - TP_m4
    Miss_m4  = length(index_true) - length(which(index_true %in% S_m4))
    FDR_m4   = FP_m4/length(S_m4)
    Power_m4 = TP_m4/length(index_true)
    
    # MSE
    theta_m4        = theta_thrshed
    theta_m4[-S_m4] = 0
    MSE_m4          = mean((theta_m4 - true_theta)^2)
  
  #------------------- Variable selection -----------------------
  # Method 5 VS-D-Lasso
    S_m5     = which(abs(theta_thrshed) > 0.1)
    TP_m5    = length(which(S_m5 %in% index_true))
    FP_m5    = length(S_m5) - TP_m5
    Miss_m5  = length(index_true) - length(which(index_true %in% S_m5))
    FDR_m5   = FP_m5/length(S_m5)
    Power_m5 = TP_m5/length(index_true)
    
    # MSE
    theta_m5        = theta_thrshed
    theta_m5[-S_m5] = 0
    MSE_m5          = mean((theta_m5 - true_theta)^2)
  
  # Method 6 VS-Lasso
    fit_lasso   = cv.ncvsurv(Phi, cbind(Y, delta), penalty = "lasso")
    theta_lasso = fit_lasso$fit$beta[, fit_lasso$min]
    
    S_m6     = which(abs(theta_lasso) > 0.1)
    TP_m6    = length(which(S_m6 %in% index_true))
    FP_m6    = length(S_m6) - length(which(index_true %in% S_m6))
    Miss_m6  = length(index_true) - length(which(index_true %in% S_m6))   
    FDR_m6   = FP_m6/length(S_m6) 
    Power_m6 = TP_m6/length(index_true)
    
    # MSE
    theta_m6        = -theta_lasso
    theta_m6[-S_m6] = 0
    MSE_m6          = mean((theta_m6 - true_theta)^2)
  
  # Method 7 VS-MCP
    fit_MCP   = cv.ncvsurv(Phi, cbind(Y,delta), penalty = "MCP")
    theta_MCP = fit_MCP$fit$beta[, fit_MCP$min]
    
    S_m7      = which(abs(theta_MCP) > 0.01)
    TP_m7     = length(which(S_m7 %in% index_true))
    FP_m7     = length(S_m7) - length(which(index_true %in% S_m7))
    Miss_m7   = length(index_true) - length(which(index_true %in% S_m7))   
    FDR_m7    = FP_m7/length(S_m7) 
    Power_m7  = TP_m7/length(index_true)
    
    # MSE
    theta_m7        = -theta_MCP
    theta_m7[-S_m7] = 0
    MSE_m7          = mean((theta_m7 - true_theta)^2)
  
  FDR_Miss_TP_FP = c(FDR_m1, FDR_m2, FDR_m3, FDR_m4, FDR_m5, FDR_m6, FDR_m7,
                     Power_m1, Power_m2, Power_m3, Power_m4, Power_m5, Power_m6, Power_m7,
                     Miss_m1, Miss_m2, Miss_m3, Miss_m4, Miss_m5, Miss_m6, Miss_m7,
                     TP_m1, TP_m2, TP_m3, TP_m4, TP_m5, TP_m6, TP_m7,
                     FP_m1, FP_m2, FP_m3, FP_m4, FP_m5, FP_m6, FP_m7,
                     MSE_propsd, MSE_m4, MSE_m5, MSE_m6, MSE_m7) 
  
  names(FDR_Miss_TP_FP) <- c("FDR_m1", "FDR_m2", "FDR_m3", "FDR_m4", "FDR_m5", "FDR_m6", "FDR_m7",
                             "Power_m1", "Power_m2", "Power_m3", "Power_m4", "Power_m5", "Power_m6", "Power_m7",
                             "Miss_m1", "Miss_m2", "Miss_m3", "Miss_m4", "Miss_m5", "Miss_m6", "Miss_m7",
                             "TP_m1", "TP_m2", "TP_m3", "TP_m4", "TP_m5", "TP_m6", "TP_m7",
                             "FP_m1", "FP_m2", "FP_m3", "FP_m4", "FP_m5", "FP_m6", "FP_m7",
                             "MSE_propsd", "MSE_m4", "MSE_m5", "MSE_m6", "MSE_m7")
  return(FDR_Miss_TP_FP)
}

LogLogistic_FDRLevel <- sapply(alpha_set, LogLogistic_fun, n = n, S = S, rho = rho, A = A,
                               d = d, q = q, censor_rate = censor_rate, Y = Y, Phi = Phi, delta = delta)

filename <- paste(seed, "LogLogistic_FDRLevel.csv",sep="")
write.csv(LogLogistic_FDRLevel %>% unlist, filename)


