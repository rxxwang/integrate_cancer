library(survival)
library(Rcpp)
library(RcppArmadillo)

#These are the same functions, that allow for multiple parameters
#But are also slower since all operations are for matrices
cppFunction(depends = "RcppArmadillo",
            ' 
  arma::mat mat_sqrt(arma::mat X) {
    arma::mat X_sqrt = sqrtmat_sympd(X); 
    return X_sqrt;
  }
'
)

fit.weighted.condlogistic <- function(weights, datasets, tol = 1e-6, trace = F) {
  p = ncol(datasets[[1]]$x)
  beta_hat <- rep(0, p)
  
  pair_terms = vector("list", length(datasets))
  
  for(i in 1:length(datasets)){
    #pair_terms[[i]]$sum = rowsum(datasets[[i]]$x, as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    temp = length(datasets[[i]]$x)
    pair_terms[[i]]$sum = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] + datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    #pair_terms[[i]]$diff = rowsum(datasets[[i]]$x * rep(c(1, -1), times = nrow(datasets[[i]]$x) / 2), 
    #                              as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    pair_terms[[i]]$diff = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] - datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    
    pair_terms[[i]]$odd_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 1)
    pair_terms[[i]]$even_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 0)
  }
  
  grad_vec = rep(0, p)
  hess = matrix(0, p, p)
  
  weights = c(1, weights)
  
  for(i in 1:length(datasets)){
    if (ncol(pair_terms[[1]]$diff) == 1 & p == 1) {
      grad_vec = grad_vec + weights[i] * sum(pair_terms[[i]]$diff * 
                                               (1 - exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) / 
                                                  (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) + 
                                                     exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta_hat))))
      hess = hess + weights[i] * sum(-pair_terms[[i]]$diff *
                                       c(exp(pair_terms[[i]]$sum %*% beta_hat) / 
                                           (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) + 
                                              exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta_hat)) ^ 2) * 
                                       pair_terms[[i]]$diff)
    } else {
      grad_vec = grad_vec + weights[i] * t(pair_terms[[i]]$diff) %*% 
        (1 - exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) / 
           (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) + 
              exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta_hat)))
      
      hess = hess + weights[i] * -t(pair_terms[[i]]$diff) %*% 
        diag(c(exp(pair_terms[[i]]$sum %*% beta_hat) / 
                 (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta_hat) + 
                    exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta_hat)) ^ 2)) %*% 
        pair_terms[[i]]$diff
    }
    
  }
  
  
  beta = beta_hat - solve(hess) %*% grad_vec
  beta_hat = cbind(beta_hat, beta)
  
  j = 2
  
  while(norm(beta_hat[, j] - beta_hat[, j - 1], type = "2") > tol){
    grad_vec = rep(0, p)
    hess = matrix(0, p, p)
    
    for(i in 1:length(datasets)){
      if (ncol(pair_terms[[1]]$diff) == 1 & p == 1) {
        grad_vec = grad_vec + weights[i] * sum(pair_terms[[i]]$diff * 
                                                 (1 - exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) / 
                                                    (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                                                       exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta))))
        hess = hess + weights[i] * sum(-pair_terms[[i]]$diff *
                                         c(exp(pair_terms[[i]]$sum %*% beta) / 
                                             (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                                                exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2) * 
                                         pair_terms[[i]]$diff)
      } else {
        grad_vec = grad_vec + weights[i] * t(pair_terms[[i]]$diff) %*% 
          (1 - exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) / 
             (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)))
        
        hess = hess + weights[i] * -t(pair_terms[[i]]$diff) %*% 
          diag(c(exp(pair_terms[[i]]$sum %*% beta) / 
                   (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                      exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)) %*% 
          pair_terms[[i]]$diff
      }
    }
    
    
    beta = beta - solve(hess) %*% grad_vec
    beta_hat = cbind(beta_hat, beta)
    
    j = j + 1
  }
  
  if(trace == T){
    return(beta_hat)
  } else{
    return(c(beta))
  }
}

test.coefs.condlogistic <- function(datasets, trace = F) {
  K = length(datasets)
  p.values = rep(1, K - 1)
  
  for (i in 1:(K - 1)) {
    strata = c(datasets[[1]]$strata, datasets[[i + 1]]$strata)
    
    n1 = length(datasets[[1]]$y)
    ni = length(datasets[[i + 1]]$y)
    p = ncol(datasets[[1]]$x)
    y1 = c(datasets[[1]]$y, datasets[[i + 1]]$y)
    x1 = rbind(datasets[[1]]$x, datasets[[i + 1]]$x)
    glm1 = clogit(y1 ~ 0 + x1 + strata(strata))
    logLik1 = logLik(glm1)
    
    y2 = c(datasets[[1]]$y, datasets[[i + 1]]$y)
    x2 = rbind(cbind(datasets[[1]]$x, matrix(0, n1, p)), cbind(matrix(0, ni, p), datasets[[i + 1]]$x))
    glm2 = clogit(y2 ~ 0 + x2 + strata(strata))
    logLik2 = logLik(glm2)
    
    LR = 2 * (logLik2 - logLik1) #deviance estimator
    p.values[i] = pchisq(LR, df = p, lower.tail = F)
  }
  return(list(p.values = p.values))
  
}

loocv.weighted.condlogistic <- function(weights, datasets, betafitted, trace = F) {
  p = ncol(datasets[[1]]$x)
  beta = betafitted
  weights = c(1, weights)
  
  pair_terms = vector("list", length(datasets))
  xtwx = matrix(0, p, p)
  
  for(i in 1:length(datasets)){
    #pair_terms[[i]]$sum = rowsum(datasets[[i]]$x, as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    temp = length(datasets[[i]]$x)
    pair_terms[[i]]$sum = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] + datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    #pair_terms[[i]]$diff = rowsum(datasets[[i]]$x * rep(c(1, -1), times = nrow(datasets[[i]]$x) / 2), 
    #                              as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    pair_terms[[i]]$diff = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] - datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    
    pair_terms[[i]]$odd_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 1)
    pair_terms[[i]]$even_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 0)
    
    if (ncol(pair_terms[[1]]$diff) == 1) {
      xtwx = xtwx + weights[i] * sum(pair_terms[[i]]$diff *
                                       c(exp(pair_terms[[i]]$sum %*% beta) / 
                                           (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                                              exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2) *
                                       pair_terms[[i]]$diff)      
    } else {
      xtwx = xtwx + weights[i] * t(pair_terms[[i]]$diff) %*%
        diag(c(exp(pair_terms[[i]]$sum %*% beta) /
                 (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) +
                    exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)) %*%
        pair_terms[[i]]$diff
    }
  }
  
  xtwx_soln = solve(xtwx)
  
  
  if (p == 1) {
    weight_loc = c(exp(pair_terms[[1]]$sum %*% beta) / 
                     (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
                        exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)
    V_matrix = as.numeric(pair_terms[[1]]$diff)^2 * weight_loc * as.numeric(xtwx_soln)
  } else {
    ##preparing local terms, because the weight matrix is diagonal, we can just use a regular square root
    ##instead of a matrix square root 
    weight_loc_sqrt = mat_sqrt(diag(c(exp(pair_terms[[1]]$sum %*% beta) / 
                                        (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
                                           exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)))
    
    V_matrix = diag(weight_loc_sqrt %*% pair_terms[[1]]$diff %*% 
                      xtwx_soln %*% t(pair_terms[[1]]$diff) %*% weight_loc_sqrt)
  }
  
  mu_local = exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) / 
    (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
       exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta))
  
  beta_loo = matrix(beta, nrow = nrow(datasets[[1]]$x) / 2, ncol = length(beta), byrow=TRUE)
  
  #vectorized way to get leave-one-out parameters
  if (ncol(pair_terms[[1]]$diff) == 1) {
    beta_loo = beta_loo - xtwx_soln[1,1] * pair_terms[[1]]$diff[,1] * (1 - mu_local)/ (1 - V_matrix)
  }else{
    beta_loo = beta_loo - t((xtwx_soln %*% t(pair_terms[[1]]$diff) %*% 
                               diag(c(1 - mu_local))))/ (1 - V_matrix)
  }  
  
  #The corresponding looped version, not used
  # for(i in 1:nrow(beta_loo)){
  #   beta_loo[i, ] = beta_loo[i, ] - c((xtwx_soln %*% t(pair_terms[[1]]$diff[i, , drop = F]) %*% 
  #                                        (1 - mu_local[i]))/ (1 - V_matrix[i]))
  # }
  
  #the cross-validated partial likelihood, taking the likelihood contribution of the strata
  #using the leave-one-out estimates
  
  #Again vectorized version of this
  if (ncol(pair_terms[[1]]$diff) == 1) {
    temp = pair_terms[[1]]$diff * beta_loo
    cvpl = sum(temp - log(1 + exp(temp)))
  } else {
    cvpl = sum(diag(pair_terms[[1]]$diff %*% t(beta_loo) - 
                      log(1 + exp(pair_terms[[1]]$diff %*% t(beta_loo)))))
  }
  
  #The corresponding looped version, not used
  # cvpl = 0
  # for(i in 1:nrow(beta_loo)){
  #   cvpl = cvpl + pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ] - 
  #     log(1 + exp(pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ]))
  # }
  #cvpl = c(cvpl)
  
  if(trace == T){
    return(list(beta_loo = beta_loo, cvpl = cvpl))
  } else{
    return(cvpl)
  }
}

loocv.weighted.condlogistic.nobeta <- function(weights, datasets, trace = F) {
  p = ncol(datasets[[1]]$x)
  beta = fit.weighted.condlogistic(weights = weights, datasets = datasets)
  
  weights = c(1, weights)
  
  pair_terms = vector("list", length(datasets))
  xtwx = matrix(0, p, p)
  
  for(i in 1:length(datasets)){
    pair_terms[[i]]$sum = rowsum(datasets[[i]]$x, as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    pair_terms[[i]]$diff = rowsum(datasets[[i]]$x * rep(c(1, -1), times = nrow(datasets[[i]]$x) / 2), 
                                  as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    
    pair_terms[[i]]$odd_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 1)
    pair_terms[[i]]$even_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 0)
    
    xtwx = xtwx + weights[i] * t(pair_terms[[i]]$diff) %*% 
      diag(c(exp(pair_terms[[i]]$sum %*% beta) / 
               (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                  exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)) %*% 
      pair_terms[[i]]$diff
  }
  
  xtwx_soln = solve(xtwx)
  
  ##preparing local terms, because the weight matrix is diagonal, we can just use a regular square root
  ##instead of a matrix square root 
  weight_loc_sqrt = mat_sqrt(diag(c(exp(pair_terms[[1]]$sum %*% beta) / 
                                      (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
                                         exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)))
  
  V_matrix = diag(weight_loc_sqrt %*% pair_terms[[1]]$diff %*% 
                    xtwx_soln %*% t(pair_terms[[1]]$diff) %*% weight_loc_sqrt)
  
  mu_local = exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) / 
    (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
       exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta))
  
  beta_loo = matrix(beta, nrow = nrow(datasets[[1]]$x) / 2, ncol = length(beta), byrow=TRUE)
  
  #vectorized way to get leave-one-out parameters
  beta_loo = beta_loo - t((xtwx_soln %*% t(pair_terms[[1]]$diff) %*% 
                             diag(c(1 - mu_local))))/ (1 - V_matrix)
  
  
  #The corresponding looped version, not used
  # for(i in 1:nrow(beta_loo)){
  #   beta_loo[i, ] = beta_loo[i, ] - c((xtwx_soln %*% t(pair_terms[[1]]$diff[i, , drop = F]) %*% 
  #                                        (1 - mu_local[i]))/ (1 - V_matrix[i]))
  # }
  
  #the cross-validated partial likelihood, taking the likelihood contribution of the strata
  #using the leave-one-out estimates
  
  #Again vectorized version of this
  cvpl = sum(diag(pair_terms[[1]]$diff %*% t(beta_loo) - 
                    log(1 + exp(pair_terms[[1]]$diff %*% t(beta_loo)))))
  
  #The corresponding looped version, not used
  # cvpl = 0
  # for(i in 1:nrow(beta_loo)){
  #   cvpl = cvpl + pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ] - 
  #     log(1 + exp(pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ]))
  # }
  #cvpl = c(cvpl)
  
  if(trace == T){
    return(list(beta_loo = beta_loo, cvpl = cvpl))
  } else{
    return(cvpl)
  }
}

loocv.weighted.condlogistic.2param = function(alphabeta, datasets, pvals, trace = F){
  p = ncol(datasets[[1]]$x)
  weights = rep(0, length(datasets) - 1)
  
  for(i in 1:length(weights)){ #getting weights from p-values
    weights[i] = ifelse(pvals[i] < alphabeta[1], 0,
                        ifelse(pvals[i] > alphabeta[2], 1, 
                               (pvals[i] - alphabeta[1]) / (alphabeta[2] - alphabeta[1])))
  }
  
  beta = fit.weighted.condlogistic(weights = weights, datasets = datasets)
  weights = c(1, weights)
  
  pair_terms = vector("list", length(datasets))
  xtwx = matrix(0, p, p)
  
  for(i in 1:length(datasets)){
    #pair_terms[[i]]$sum = rowsum(datasets[[i]]$x, as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    temp = length(datasets[[i]]$x)
    pair_terms[[i]]$sum = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] + datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    #pair_terms[[i]]$diff = rowsum(datasets[[i]]$x * rep(c(1, -1), times = nrow(datasets[[i]]$x) / 2), 
    #                              as.integer(gl(nrow(datasets[[i]]$x), 2, nrow(datasets[[i]]$x))))
    pair_terms[[i]]$diff = as.matrix(datasets[[i]]$x[seq(1, temp - 1, by = 2)] - datasets[[i]]$x[seq(2, temp, by = 2)], ncol = 1)
    
    pair_terms[[i]]$odd_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 1)
    pair_terms[[i]]$even_indices = which(1:nrow(datasets[[i]]$x) %% 2 == 0)
    
    if (ncol(pair_terms[[1]]$diff) == 1) {
      xtwx = xtwx + weights[i] * sum(pair_terms[[i]]$diff *
                                       c(exp(pair_terms[[i]]$sum %*% beta) / 
                                           (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) + 
                                              exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2) *
                                       pair_terms[[i]]$diff)      
    } else {
      xtwx = xtwx + weights[i] * t(pair_terms[[i]]$diff) %*%
        diag(c(exp(pair_terms[[i]]$sum %*% beta) /
                 (exp(datasets[[i]]$x[pair_terms[[i]]$odd_indices, , drop = FALSE] %*% beta) +
                    exp(datasets[[i]]$x[pair_terms[[i]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)) %*%
        pair_terms[[i]]$diff
    }
  }
  
  xtwx_soln = solve(xtwx)
  
  if (p == 1) {
    weight_loc = c(exp(pair_terms[[1]]$sum %*% beta) / 
                     (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
                        exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)
    V_matrix = as.numeric(pair_terms[[1]]$diff)^2 * weight_loc * as.numeric(xtwx_soln)
  } else {
    ##preparing local terms, because the weight matrix is diagonal, we can just use a regular square root
    ##instead of a matrix square root 
    weight_loc_sqrt = mat_sqrt(diag(c(exp(pair_terms[[1]]$sum %*% beta) / 
                                        (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
                                           exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta)) ^ 2)))
    
    V_matrix = diag(weight_loc_sqrt %*% pair_terms[[1]]$diff %*% 
                      xtwx_soln %*% t(pair_terms[[1]]$diff) %*% weight_loc_sqrt)
  }
  
  mu_local = exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) / 
    (exp(datasets[[1]]$x[pair_terms[[1]]$odd_indices, , drop = FALSE] %*% beta) + 
       exp(datasets[[1]]$x[pair_terms[[1]]$even_indices, , drop = FALSE] %*% beta))
  
  beta_loo = matrix(beta, nrow = nrow(datasets[[1]]$x) / 2, ncol = length(beta), byrow=TRUE)
  
  #vectorized way to get leave-one-out parameters
  if (ncol(pair_terms[[1]]$diff) == 1) {
    beta_loo = beta_loo - xtwx_soln[1,1] * pair_terms[[1]]$diff[,1] * (1 - mu_local)/ (1 - V_matrix)
  }else{
    beta_loo = beta_loo - t((xtwx_soln %*% t(pair_terms[[1]]$diff) %*% 
                               diag(c(1 - mu_local))))/ (1 - V_matrix)
  }
  
  
  #The corresponding looped version, not used
  # for(i in 1:nrow(beta_loo)){
  #   beta_loo[i, ] = beta_loo[i, ] - c((xtwx_soln %*% t(pair_terms[[1]]$diff[i, , drop = F]) %*% 
  #                                        (1 - mu_local[i]))/ (1 - V_matrix[i]))
  # }
  
  #the cross-validated partial likelihood, taking the likelihood contribution of the strata
  #using the leave-one-out estimates
  
  #Again vectorized version of this
  if (ncol(pair_terms[[1]]$diff) == 1) {
    temp = pair_terms[[1]]$diff * beta_loo
    cvpl = sum(temp - log(1 + exp(temp)))
  } else {
    cvpl = sum(diag(pair_terms[[1]]$diff %*% t(beta_loo) - 
                      log(1 + exp(pair_terms[[1]]$diff %*% t(beta_loo)))))
  
  #The corresponding looped version, not used
  # cvpl = 0
  # for(i in 1:nrow(beta_loo)){
  #   cvpl = cvpl + pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ] - 
  #     log(1 + exp(pair_terms[[1]]$diff[i, ] %*% beta_loo[i, ]))
  # }
  #cvpl = c(cvpl)
  
  if(trace == T){
    return(list(beta_loo = beta_loo, cvpl = cvpl))
  } else{
    return(cvpl)
  }
  }
}

joint.fitting = function(datasets, fit.weighted = fit.weighted.condlogistic, 
                         loocv.weighted = loocv.weighted.condlogistic, 
                         test.coefs = test.coefs.condlogistic, 
                         method = "testingopt", trace = F) {
  #print(Sys.time())
  K = length(datasets)
  best.weights = rep(0, K-1) #initialize weights as a vector of 0s
  
  if (method == "opt"){ #full optimization, null hypothesis of no weights
    
    opt = optim(par = rep(0.5, K - 1), fn = loocv.weighted.condlogistic.nobeta, datasets = datasets, method = "L-BFGS-B", 
                lower = rep(0, K - 1), upper = rep(1, K - 1), control = list(fnscale = -1))
    best.weights = opt$par
    loocv = opt$value
    
  } else if (method == "testingopt"){ #optimize the 2 param loo, minimizing the loocv value
    ##use the better of binary or pval for the starting point
    pvals = as.numeric(test.coefs(datasets)$p.values)
    loocvpval = loocv.weighted(pvals, datasets, fit.weighted(pvals, datasets))
    pvals.bin = as.numeric(pvals > 0.05)
    loocvbin = loocv.weighted(pvals.bin, datasets, fit.weighted(pvals.bin, datasets))
    
    if(loocvpval > loocvbin){ #pval is better
      opt = optim(par = c(0, 1), fn = loocv.weighted.condlogistic.2param, datasets = datasets, pvals = pvals,
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1), control = list(fnscale = -1, factr = 1e9, pgtol = 1e-3, maxit = 10))
      #opt = optim(par = c(0, 1), fn = loocv.weighted.condlogistic.2param, datasets = datasets, pvals = pvals,
      #                   method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1), control = list(fnscale = -1))
    } else {
      opt = optim(par = c(0.05, 0.05), fn = loocv.weighted.condlogistic.2param, datasets = datasets, pvals = pvals,
                         method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1), control = list(fnscale = -1, factr = 1e9, pgtol = 1e-3, maxit = 10))
      #opt = optim(par = c(0.05, 0.05), fn = loocv.weighted.condlogistic.2param, datasets = datasets, pvals = pvals,
      #                  method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1), control = list(fnscale = -1))
    }
    testingopt = opt$par
    #print(opt)
    
    for(i in 1:(length(best.weights))){#convert p-values to weights
      
      best.weights[i] = ifelse(pvals[i] < testingopt[1], 0,
                               ifelse(pvals[i] > testingopt[2], 1, (pvals[i] - testingopt[1])/
                                        (testingopt[2] - testingopt[1])))
      best.weights[i] = ifelse(sum(datasets[[i+1]]$x) == 0, 0, best.weights[i])

    }
    
    loocv = loocv.weighted(best.weights, datasets, fit.weighted(best.weights, datasets))
  }
  
  if (trace) {
    print(method)
    print("best.weights")
    print(best.weights)
  }
  
  if(method == "opt" | method == "testingopt"){
    
    return (list(weights = best.weights, par = fit.weighted(best.weights, datasets), 
                 #par.localonly = fit.weighted(rep(0, K-1), datasets), 
                 opt = opt,
                 loocv = loocv))
  }
}