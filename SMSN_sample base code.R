library(sn)
library(moments)
library(histogram)
library(nnls)
library(tidyverse)
library(nspmix)
library(plot3D)
library(mvtnorm)
library(pbivnorm)
# MVT and MST
library(fitHeavyTail)
# library(EMMIXuskew)
# ARI AMI
library(aricode)

pbivnorm_pmvnormlike = function(lower, mu, Sigma)
{
  if (length(mu) == 1)
  {
    result = 1 - pnorm(q = lower, mean = mu, sd = sqrt(Sigma))
    return(result)
  }
  
  Sigma_diag_sqrt = (1/(Sigma %>% diag() %>% sqrt())) %>% diag()
  
  lower_std = (lower - mu) %*% Sigma_diag_sqrt
  upper_std = -lower_std
  
  Corr_mat = Sigma_diag_sqrt %*% Sigma %*% Sigma_diag_sqrt
  rho = Corr_mat[1, 2]
  
  result = pbivnorm(x = upper_std, rho = rho)
  
  # upper_std_temp = upper_std
  # iter_result = 0
  # while ((is.na(result)) & (iter_result < 100))
  # {
  #   iter_result = iter_result + 1
  #   upper_std_temp = upper_std_temp*(0.9)
  #   result = pbivnorm(x = upper_std_temp, rho = rho)
  # }
  
  return(result)
}

dmsn_sahu = function(x, xi, Sigma, Lambda) #Lin, 2009
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  
  if (!is.null(dim(x)))
  {
    result = 
      apply(x, 1, FUN = dmsn_sahu, 
            xi = xi, Sigma = Sigma, Lambda = Lambda)
  }else{
    p = length(x)
    
    #Lambda = diag(lambda)
    Omega = Sigma + (Lambda %*% t(Lambda))
    Omega_inv = solve(Omega)
    Delta = diag(p) - (t(Lambda) %*% Omega_inv %*% Lambda)
    
    constant_p = 2^p
    phi_p = dmvnorm(x = x, mean = xi, sigma = Omega)
    Phi_p = pmvnorm(upper = as.numeric(t(Lambda) %*% Omega_inv %*% (x - xi)), sigma = Delta) %>% as.numeric()
    
    result = constant_p * phi_p * Phi_p 
  }
  
  return(result)
}

dmsn_sahu_v2 = function(x, xi, Sigma, Lambda, 
                        Omega = Sigma + (Lambda %*% t(Lambda)),
                        Omega_inv = solve(Omega)) #Lin, 2009
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  
  p = length(xi)
  
  Delta = diag(p) - (t(Lambda) %*% Omega_inv %*% Lambda)
  
  xi_mat = matrix(data = xi, nrow = nrow(x), ncol = p, byrow = T)
  upper_mat = t(Lambda) %*% Omega_inv %*% t(x - xi_mat) %>% t()
  
  constant_p = 2^p
  phi_p = dmvnorm(x = x, mean = xi, sigma = Omega)
  Phi_p = 
    apply(X = upper_mat,
          MARGIN = 1,
          FUN = function(upper_vec){pmvnorm(upper = upper_vec, sigma = Delta)}
    )
  
  result = constant_p * phi_p * Phi_p 
  
  return(result)
}

dmsn_sahu_v3_biv = function(x, xi, Sigma, Lambda, 
                            Omega = Sigma + (Lambda %*% t(Lambda)),
                            Omega_inv = solve(Omega)) #Lin, 2009
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  
  p = length(xi)
  
  Delta = diag(p) - (t(Lambda) %*% Omega_inv %*% Lambda)
  
  xi_mat = matrix(data = xi, nrow = nrow(x), ncol = p, byrow = T)
  upper_mat = t(Lambda) %*% Omega_inv %*% t(x - xi_mat) %>% t()
  
  constant_p = 2^p
  phi_p = dmvnorm(x = x, mean = xi, sigma = Omega)
  # Phi_p =
  #   apply(X = upper_mat,
  #         MARGIN = 1,
  #         FUN = function(x){pmvnorm(upper = x, sigma = Delta, algorithm = Miwa())}
  #   )
  Delta_diag_sqrt = sqrt(diag(Delta))
  Delta_diag_outer = outer(X = Delta_diag_sqrt, Y = Delta_diag_sqrt)
  Delta_rho = Delta/Delta_diag_outer
  Delta_diag_sqrt_mat = matrix(Delta_diag_sqrt, 
                               nrow = nrow(upper_mat),
                               ncol = ncol(upper_mat), byrow = T)
  upper_mat_normalized = upper_mat/Delta_diag_sqrt_mat
  Phi_p = pbivnorm(x = upper_mat_normalized, rho = Delta_rho[1,2])
  
  result = constant_p * phi_p * Phi_p 
  result[(result > -1e-16) & (result < 0)] = 0
  
  return(result)
}

dmst_sahu = function(x, xi, Sigma, Lambda, nu,
                     Omega = Sigma + (Lambda %*% t(Lambda)),
                     Omega_inv = solve(Omega))
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  
  n = nrow(x)
  p = nrow(Sigma)
  
  Delta = diag(p) - (t(Lambda) %*% Omega_inv %*% Lambda)
  xi_mat = matrix(xi, nrow = n, ncol = p, byrow = T)
  
  dx = 
    (x - xi_mat) %>%
    apply(MARGIN = 1, FUN = function(x_vec){t(x_vec) %*% Omega_inv %*% x_vec})
  
  q = t(t(Lambda) %*% Omega_inv %*% t(x - xi_mat))
  q_2 = sqrt((nu + p)/(nu + dx)) %>% matrix(nrow = n, ncol = p)
  
  x_star = q*q_2
  
  constant_p = 2^p
  t_p = dmvt(x = x, delta = xi, sigma = Omega, df = nu, log = F)
  T_p =
    x_star %>%
    apply(MARGIN = 1,
          FUN = function(upper_vec){pmvt(upper = upper_vec,
                                         delta = rep(0, times = p), sigma = Delta, df = round(nu))})
  
  result = constant_p*t_p*T_p
  
  return(result)
}

rmsn_sahu = function(n, xi, Sigma, Lambda) #Lin, 2009
{
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  
  p = nrow(Sigma)
  
  tau = rmvnorm(n = n, sigma = diag(p)) %>% abs() #half normal
  U = rmvnorm(n = n, sigma = Sigma)
  xi_mat = matrix(xi, nrow = n, ncol = length(xi), byrow = T)
  
  result = xi_mat + (tau %*% t(Lambda)) + U
  
  return(result)
}

rmst_sahu = function(n, xi, Sigma, Lambda, nu) #Lee, 2011
{
  if ( is.null(dim(Lambda)) )
  {
    Lambda = diag(Lambda)
  } else {
    Lambda = Lambda
  }
  p = nrow(Sigma)
  
  I_mat = diag(p)
  
  w = rgamma(n = n, shape = nu/2, rate = nu/2)
  w_mat = matrix(w, nrow = n, ncol = p)
  U = abs(rmvnorm(n = n, mean = rep(0, times = p), sigma = I_mat)/sqrt(w_mat))
  xi_mat = matrix(xi, nrow = n, ncol = p, byrow = T)
  
  Y = (rmvnorm(n = n, sigma = Sigma)/sqrt(w_mat)) + (U %*% t(Lambda) + xi_mat)
  
  return(Y)
}

rsmsn_sahu = function(n, Q_support)
{
  n_set = round(n*Q_support$prob)
  target_idx = sample(x = 1:length(n_set), size = 1)
  
  n_set[target_idx] = n - sum(n_set[-target_idx])
  
  n_alpha_mat = cbind(n_set, Q_support$alpha)
  
  colnames(n_alpha_mat) = c("n_temp", "alpha_temp")
  rownames(n_alpha_mat) = str_c("sample", 1:length(n_set))
  
  n_alpha_list = n_alpha_mat %>% t() %>% as.data.frame() %>% as.list()
  
  rsmsn_sample =
    lapply(X = n_alpha_list,
          FUN = function(n_alpha_vec){
            rmsn_sahu(n = n_alpha_vec[1],
                      xi = Q_support$xi, 
                      Sigma = n_alpha_vec[2]*Q_support$Sigma,
                      Lambda = Q_support$Lambda)
            
            }) 
  
  rsmsn_sample = do.call(what = rbind, args = rsmsn_sample)
  
  return(rsmsn_sample)
}

# 2. Calculate density of mixture distribution, Q_support is a list containing parameters

density_mixture = function(x, Q_support)
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  num_support = length(Q_support$alpha)
  
  result = 0
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  for (i in 1:num_support)
  {
    result = result +
      Q_support$prob[i]*dmsn_sahu_v3_biv(x,
                                          xi = Q_xi,
                                          Sigma = Q_support$alpha[i]*Q_Sigma,
                                          Lambda = Q_Lambda)
  }
  
  return(result)
}

density_mixture_v2 = function(x, Q_support)
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  num_support = length(Q_support$alpha)
  
  result = 0
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  Q_prob_alpha_mat = cbind(Q_support$prob, Q_support$alpha)
  
  Q_prob_alpha_list = 
    split(Q_prob_alpha_mat, 
          rep(1:nrow(Q_prob_alpha_mat), times = ncol(Q_prob_alpha_mat)))
  
  result_list =
    lapply(Q_prob_alpha_list, 
           FUN = function(prob_alpha) #corrected
           {
             result = prob_alpha[1]*dmsn_sahu_v3_biv(x = x,
                                                      xi = Q_xi,
                                                      Sigma = prob_alpha[2]*Q_Sigma,
                                                      Lambda = Q_Lambda)
             return(result)
           })
  n = nrow(x)
  
  result_mat = matrix(unlist(result_list), ncol = n, byrow = T)
  result = apply(result_mat, 2, FUN = sum)
  
  return(result)
}

density_mixture_v2_mst_sahu = function(x, Q_support) 
{
  if (is.null(dim(x))) {x = matrix(x, nrow = 1)}
  
  num_support = length(Q_support$alpha)
  
  result = 0
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  Q_nu = Q_support$nu
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  Q_prob_alpha_mat = cbind(Q_support$prob, Q_support$alpha)
  
  Q_prob_alpha_list = 
    split(Q_prob_alpha_mat, 
          rep(1:nrow(Q_prob_alpha_mat), times = ncol(Q_prob_alpha_mat)))
  
  result_list =
    lapply(Q_prob_alpha_list, #corrected
           FUN = function(prob_alpha)
           {
             result = prob_alpha[1]*dmst_sahu(x = x,
                                              xi = Q_xi,
                                              Sigma = prob_alpha[2]*Q_Sigma,
                                              Lambda = Q_Lambda,
                                              nu = Q_nu)
             return(result)
           })
  n = nrow(x)
  
  result_mat = matrix(unlist(result_list), ncol = n, byrow = T)
  result = apply(result_mat, 2, FUN = sum)
  
  return(result)
}

# 3. Directional Derivative function and log likelihood function which is used to calculate optimal mass that new support should have

Directional = function(component_parameter, Q_support, x)
{
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  denominator_likelihood = density_mixture_v2(x, Q_support = Q_support)
  
  result = ( (dmsn_sahu_v3_biv(x, 
                        xi = Q_xi,
                        Sigma = component_parameter*Q_Sigma,
                        Lambda = Q_Lambda) / denominator_likelihood) %>% sum() )- n
  
  return(result)
}

Directional_with_identifiability = function(component_parameter, Q_support, x, idx_support_to_be_kept = 1)
{
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  Q_support_to_be_kept = Q_support
  Q_support_not_to_be_kept = Q_support
  
  Q_support_to_be_kept$alpha = Q_support$alpha[idx_support_to_be_kept]
  Q_support_to_be_kept$prob = Q_support$prob[idx_support_to_be_kept]
  
  idx_support_not_to_be_kept = ((1:length(Q_support$alpha)) != idx_support_to_be_kept)
  
  Q_support_not_to_be_kept$alpha = Q_support$alpha[idx_support_not_to_be_kept]
  Q_support_not_to_be_kept$prob = Q_support$prob[idx_support_not_to_be_kept]
  
  if (sum(idx_support_not_to_be_kept) == 0)
  {
    Q_support_not_to_be_kept$alpha = Q_support_to_be_kept$alpha
    Q_support_not_to_be_kept$prob = 0
  }
  
  pi_tbk = Q_support$prob[idx_support_to_be_kept]
  
  denominator_likelihood_1 = 
    (pi_tbk *
     dmsn_sahu_v3_biv(x = x, 
                      xi = Q_xi,
                      Sigma = Q_support_to_be_kept$alpha * Q_Sigma,
                      Lambda = Q_Lambda))
  denominator_likelihood_2 = density_mixture_v2(x, Q_support = Q_support_not_to_be_kept)
  
  numerator_likelihood_1 = 
    (pi_tbk * 
       dmsn_sahu_v3_biv(x, 
                 xi = Q_xi, 
                 Sigma = Q_support_to_be_kept$alpha*Q_Sigma,
                 Lambda = Q_Lambda))
  
  numerator_likelihood_component = 
    ((1 - pi_tbk) * 
       dmsn_sahu_v3_biv(x,
                 xi = Q_xi,
                 Sigma = component_parameter*Q_Sigma,
                 Lambda = Q_Lambda))
  
  numerator_likelihood = numerator_likelihood_1 + numerator_likelihood_component
  denominator_likelihood = denominator_likelihood_1 + denominator_likelihood_2
  
  result = ( (numerator_likelihood / denominator_likelihood) %>% sum() )- n
  
  return(result)
}

L_function_vem = function(mixing_p, maximizer_parameter, minimizer_idx, Q_support, x)
{
  P_likelihood = density_mixture_v2(x, Q_support = Q_support)
  
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  result = P_likelihood
  result = result - 
    (mixing_p)*Q_support$prob[minimizer_idx]*dmsn_sahu_v3_biv(x,
                                                       xi = Q_xi,
                                                       Sigma = Q_support$alpha[minimizer_idx]*Q_Sigma,
                                                       Lambda = Q_Lambda) +
    (mixing_p)*Q_support$prob[minimizer_idx]*dmsn_sahu_v3_biv(x,
                                                       xi = Q_xi,
                                                       Sigma = maximizer_parameter*Q_Sigma,
                                                       Lambda = Q_Lambda)
  
  result = sum(log(result))
  
  return(result)
}

L_function_vem_with_identifiability = function(mixing_beta, maximizer_parameter, minimizer_idx, Q_support, x, idx_support_to_be_kept = 1)
{
  Q_xi = Q_support$xi
  Q_Sigma = Q_support$Sigma
  
  if ( is.null(dim(Q_support$Lambda)) )
  {
    Q_Lambda = diag(Q_support$Lambda)
  } else {
    Q_Lambda = Q_support$Lambda
  }
  
  Q_support_to_be_kept = Q_support
  Q_support_not_to_be_kept = Q_support
  
  Q_support_to_be_kept$alpha = Q_support$alpha[idx_support_to_be_kept]
  Q_support_to_be_kept$prob = Q_support$prob[idx_support_to_be_kept]
  
  idx_support_not_to_be_kept = ((1:length(Q_support$alpha)) != idx_support_to_be_kept)
  
  Q_support_not_to_be_kept$alpha = Q_support$alpha[idx_support_not_to_be_kept]
  Q_support_not_to_be_kept$prob = Q_support$prob[idx_support_not_to_be_kept]
  
  if (sum(idx_support_not_to_be_kept) == 0)
  {
    Q_support_not_to_be_kept$alpha = Q_support_to_be_kept$alpha
    Q_support_not_to_be_kept$prob = 0
  }
  
  pi_tbk = Q_support$prob[idx_support_to_be_kept]
  
  result_1 = pi_tbk * dmsn_sahu_v3_biv(x = x,
                                xi = Q_xi,
                                Sigma = Q_support_to_be_kept$alpha * Q_Sigma,
                                Lambda = Q_Lambda)
  
  result_2 = density_mixture_v2(x = x, Q_support = Q_support_not_to_be_kept) #(1 - pi_tbk) are already in Q_support_not_to_be_kept
  
  result_3 = -(1 - pi_tbk) * mixing_beta * Q_support$prob[minimizer_idx] * dmsn_sahu_v3_biv(x = x, 
                                                                                     xi = Q_xi,
                                                                                     Sigma = Q_support$alpha[minimizer_idx]*Q_Sigma,
                                                                                     Lambda = Q_Lambda)
  
  result_4 = (1 - pi_tbk) * mixing_beta * Q_support$prob[minimizer_idx] * dmsn_sahu_v3_biv(x = x,
                                                                                    xi = Q_xi,
                                                                                    Sigma = maximizer_parameter*Q_Sigma,
                                                                                    Lambda = Q_Lambda)
  
  result = result_1 + result_2 + result_3 + result_4
  
  result = sum(log(result))
  
  return(result)
}

# 4. Functions for EM algorithm in Lin 2009

Theta_list_generator = function(Q_support, Z = NULL, eta_ij = NULL, Psi_ij = NULL)
{
  p = length(Q_support$xi)
  
  num_support = length(Q_support$alpha)
  
  Theta_list = vector("list", num_support)
  for (i in 1:num_support)
  {
    Theta_list[[i]]$xi = Q_support$xi
    Theta_list[[i]]$Sigma = Q_support$Sigma * Q_support$alpha[i]
    Theta_list[[i]]$Lambda = Q_support$Lambda
    Theta_list[[i]]$prob = Q_support$prob[i]
    if (!is.null(Z)) { Theta_list[[i]]$Z_i = Z[, i] }
    if (!is.null(eta_ij)) { Theta_list[[i]]$eta_i = eta_ij[[i]] }
    if (!is.null(Psi_ij)) 
    {
      Theta_list[[i]]$Psi_i =
        Psi_ij[[i]] %>% 
        unlist() %>%
        array(dim = c(p, p, length(Psi_ij[[i]])))
    }
    Theta_list[[i]]$alpha_i = Q_support$alpha[i]
  }
  
  return(Theta_list)
}

alpha_TN = function(a_vector, mu, Sigma)
{
  #result = pmvnorm(lower = a_vector, mean = mu, sigma = Sigma) %>% as.numeric()
  result = pbivnorm_pmvnormlike(lower = a_vector, mu = mu, Sigma = Sigma)
  return(result)
}

f_r = function(a_vector, mu, Sigma, r)
{
  if ( is.null(dim(a_vector)) )
  {
    a_r = a_vector[r]
  } else {
    a_r = a_vector[r, 1]
  }
  
  mu_r = mu[r]
  sigma_rr = Sigma[r,r]
  
  result = dnorm(x = a_r, mean = mu_r, sd = sqrt(sigma_rr))
  return(result)               
}

G_r = function(a_vector, mu, Sigma, r)
{
  if ( is.null(dim(a_vector)) )
  {
    a_vector_r = matrix(a_vector[-r], ncol = 1)
    a_r = matrix(a_vector[r], ncol = 1)
  } else {
    a_vector_r = a_vector[-r, 1, drop = F]
    a_r = a_vector[r, 1, drop = F]
  }
  
  mu_1 = matrix(mu[-r], ncol = 1)
  mu_2 = matrix(mu[r], ncol = 1)
  
  Sigma_11 = Sigma[-r, -r, drop = F]
  Sigma_12 = Sigma[-r, r, drop = F]
  Sigma_21 = Sigma[r, -r, drop = F]
  Sigma_22 = Sigma[r, r, drop = F]
  
  mu_c = mu_1 + ( Sigma_12 %*% (1/Sigma_22) %*% (a_r - mu_2) )
  Sigma_c = Sigma_11 - ( Sigma_12 %*% (1/Sigma_22) %*% Sigma_21 )
  
  
  result = pbivnorm_pmvnormlike(lower = as.vector(a_vector_r), mu = as.vector(mu_c), Sigma = Sigma_c)
  
  return(result)
}

f_rs = function(a_vector, mu, Sigma, r, s)
{
  if ( is.null(dim(a_vector)) )
  {
    a_r = a_vector[r]; a_s = a_vector[s]
  } else {
    a_r = a_vector[r, 1]; a_s = a_vector[s, 1]
  }
  mu_rs = mu[c(r,s)]
  Sigma_rs = Sigma[c(r,s), c(r,s)]
  
  result = dmvnorm(x = c(a_r, a_s), mean = mu_rs, sigma = Sigma_rs) %>% as.numeric()
  return(result)
}

G_rs = function(a_vector, mu, Sigma, r, s)
{
  rs = c(r, s)
  
  if ( is.null(dim(a_vector)) )
  {
    a_vector_rs = matrix(a_vector[-rs], ncol = 1)
    a_rs = matrix(a_vector[rs], ncol = 1)
  } else {
    a_vector_rs = a_vector[-rs, 1, drop = F]
    a_rs = a_vector[rs, 1]
  }
  
  mu_1 = matrix(mu[-rs], ncol = 1)
  mu_2 = matrix(mu[rs], ncol = 1)
  
  Sigma_11 = Sigma[-rs, -rs, drop = F]
  Sigma_12 = Sigma[-rs, rs, drop = F]
  Sigma_21 = Sigma[rs, -rs, drop = F]
  Sigma_22 = Sigma[rs, rs, drop = F]
  
  mu_c = mu_1 + ( Sigma_12 %*% solve(Sigma_22) %*% (a_rs - mu_2) )
  Sigma_c = Sigma_11 - ( Sigma_12 %*% solve(Sigma_22) %*% Sigma_21 )
  
  if (length(rs) == 2)
  {
    result = 1
  } else {
    
    result = pbivnorm_pmvnormlike(lower = as.vector(a_vector_r), mean = as.vector(mu_c), sigma = Sigma_c)
  }
  return(result)
}

EX_TN = function(a_vector, mu, Sigma, tol_alpha_TN = 1e-16)
{
  alpha_TN_value = as.numeric(alpha_TN(a_vector = a_vector, mu = mu, Sigma = Sigma))
  
  r_set = 1:length(a_vector)
  
  q_vector = 
    sapply(X = r_set, FUN = function(r) f_r(a_vector = a_vector, mu = mu, Sigma = Sigma, r = r)*
             G_r(a_vector = a_vector, mu = mu, Sigma = Sigma, r = r)) %>%
    matrix(ncol = 1)
  
  mu_vector = matrix(mu, ncol = 1)
  
  if (alpha_TN_value < tol_alpha_TN)
  {
    result = mu_vector
  } else {
    result = mu_vector + (1/alpha_TN_value) * (Sigma %*% q_vector)
    result = as.numeric(result)
  }
  
  return(result)
}

EXXt_TN = function(a_vector, mu, Sigma, tol_alpha_TN = 1e-16)
{
  alpha_TN_value = as.numeric(alpha_TN(a_vector = a_vector, mu = mu, Sigma = Sigma))
  eta = EX_TN(a_vector = a_vector, mu = mu, Sigma = Sigma) %>% matrix(ncol = 1)
  
  # Making matrix H
  r_set = 1:length(a_vector)
  i_mat = matrix(r_set, 
                 nrow = length(a_vector),
                 ncol = length(a_vector))
  j_mat = t(i_mat)
  i_idx = i_mat[upper.tri(i_mat)]
  j_idx = j_mat[upper.tri(j_mat)]
  ij_idx = cbind(i_idx, j_idx)
  
  H_element = 
    apply(ij_idx, 1,
          FUN = function(rs) 
          {
            f_rs(a_vector = a_vector, mu = mu, Sigma = Sigma, r = rs['i_idx'], s = rs['j_idx'])*
              G_rs(a_vector = a_vector, mu = mu, Sigma = Sigma, r = rs['i_idx'], s = rs['j_idx'])
          }
    )
  H_u = matrix(0, 
               nrow = length(a_vector),
               ncol = length(a_vector))
  H_u[upper.tri(H_u)] = H_element
  
  H = H_u + t(H_u)
  
  # Making matrix D
  Sigma_diag_inv = 1/diag(Sigma)
  Sigma_H_diag = diag(Sigma %*% H)
  
  r_set = 1:length(a_vector)
  q_vector = #f_ar * G_ar
    sapply(X = r_set, FUN = function(r) f_r(a_vector = a_vector, mu = mu, Sigma = Sigma, r = r)*
             G_r(a_vector = a_vector, mu = mu, Sigma = Sigma, r = r))
  
  D = ( Sigma_diag_inv * ( ((a_vector - mu) * q_vector) - Sigma_H_diag) ) %>% diag()
  
  #Calculate EXXt and return
  mu_erec = matrix(mu, ncol = 1)
  
  if (alpha_TN_value < tol_alpha_TN)
  {
    result = (mu_erec %*% t(eta)) + (eta %*% t(mu_erec)) - (mu_erec %*% t(mu_erec)) + Sigma
  } else {
    result = 
      (mu_erec %*% t(eta)) + (eta %*% t(mu_erec)) - (mu_erec %*% t(mu_erec)) + 
      Sigma + ( (1/alpha_TN_value) * (Sigma %*% (H + D) %*% Sigma) )
  }
  
  return(result)
}

mu_tau_TN_generator = function(Theta, x)
{
  temp_Theta = Theta
  
  dim_length = length(temp_Theta$xi)
  I_p = diag(1, ncol = dim_length, nrow = dim_length)
  
  Omega_i = temp_Theta$Sigma + ( temp_Theta$Lambda %*% t(temp_Theta$Lambda) )
  Omega_i_inv = solve(Omega_i)
  Delta_i = I_p - ( t(temp_Theta$Lambda) %*% Omega_i_inv %*% temp_Theta$Lambda )
  
  y = t(x)
  xi_mat = matrix(temp_Theta$xi, nrow = nrow(y), ncol = ncol(y))
  
  mu_tau_TN_i = ( t(temp_Theta$Lambda) %*% Omega_i_inv %*% (y - xi_mat) ) %>% t()
  
  mu_Sigma_tau_TN_i = list(mu = mu_tau_TN_i, Sigma = Delta_i)
  
  return(mu_Sigma_tau_TN_i)  
}

eta_generator = function(mu_tau_TN_i)
{
  dim_length = ncol(mu_tau_TN_i$mu)
  a_vector = rep(0, times = dim_length)
  
  mu_tau_TN_i_mu_list = split(mu_tau_TN_i$mu, rep(1:nrow(mu_tau_TN_i$mu), times = ncol(mu_tau_TN_i$mu)))
  
  eta_i = 
    sapply(X = mu_tau_TN_i_mu_list, 
           FUN = function(mu_j) EX_TN(a_vector = a_vector, mu = mu_j, Sigma = mu_tau_TN_i$Sigma)) %>% t()
  
  return(eta_i)
}

Psi_generator = function(mu_tau_TN_i)
{
  dim_length = ncol(mu_tau_TN_i$mu)
  a_vector = rep(0, times = dim_length)
  
  mu_tau_TN_i_mu_list = split(mu_tau_TN_i$mu, rep(1:nrow(mu_tau_TN_i$mu), times = ncol(mu_tau_TN_i$mu)))
  
  Psi_i = 
    lapply(X = mu_tau_TN_i_mu_list, 
           FUN = function(mu_j) EXXt_TN(a_vector = a_vector, mu = mu_j, Sigma = mu_tau_TN_i$Sigma))
  
  return(Psi_i)
}

Z_generator = function(x, Theta_list)
{
  num_support = length(Theta_list)
  
  Z_numerator = 
    lapply(X = Theta_list,
           FUN = function(Theta_i)
           {
             Theta_i$prob*dmsn_sahu_v3_biv(x = x,
                                           xi = Theta_i$xi,
                                           Sigma = Theta_i$Sigma,
                                           Lambda = Theta_i$Lambda)
           }) %>%
    unlist %>%
    matrix(nrow = nrow(x), ncol = num_support)
  
  Z_denominator = apply(Z_numerator, 1, sum) %>% matrix(ncol = ncol(Z_numerator),
                                                        nrow = nrow(Z_numerator))
  
  Z = Z_numerator / Z_denominator
  
  return(Z)
}

xi_generator = function(Theta_list, x, Z, eta_ij, Q_support)
{
  num_support = length(Q_support$alpha)
  
  numerator_term = 0
  denominator_term = 0
  
  for (i in 1:num_support)
  {
    Theta_i = Theta_list[[i]]
    Z_i_div_alpha = Z[, i] / Q_support$alpha[i]
    
    eta_i = eta_ij[[i]]
    Lambda_i = Theta_i$Lambda
    
    vector_term = t(x) - (Lambda_i %*% t(eta_i))
    scalar_term_mat = matrix(Z_i_div_alpha, 
                             ncol = ncol(vector_term), nrow = nrow(vector_term), byrow = T)
    
    numerator_i = (scalar_term_mat * vector_term) %>% apply(MARGIN = 1, FUN = sum)
    denominator_i = sum(Z_i_div_alpha)
    
    numerator_term = numerator_term + numerator_i
    denominator_term = denominator_term + denominator_i
  }
  
  result = numerator_term/denominator_term
  return(result)
}

xi_generator_v2 = function(Theta_list, x)
{
  num_support = length(Theta_list)
  
  num_denom_list =
    lapply(X = Theta_list, 
           FUN = 
            function(Theta_i) 
            {
              Lambda_i = Theta_i$Lambda
              
              Z_i_div_alpha = Theta_i$Z_i / Theta_i$alpha_i
              eta_i = Theta_i$eta_i
              
              vector_term = t(x) - (Lambda_i %*% t(eta_i))
              scalar_term_mat = matrix(Z_i_div_alpha, 
                                       ncol = ncol(vector_term), nrow = nrow(vector_term), byrow = T)
              
              numerator_i = (scalar_term_mat * vector_term) %>% apply(MARGIN = 1, FUN = sum)
              denominator_i = sum(Z_i_div_alpha)
              
              result = c(numerator_i, denominator_i)
              return(result)
            })
  
  num_denom_mat = 
    unlist(num_denom_list) %>% 
    matrix(byrow = T, nrow = num_support)
  
  num_denom_sum = apply(num_denom_mat, 2, sum)
  
  result = num_denom_sum[-length(num_denom_sum)]/num_denom_sum[length(num_denom_sum)]
  
  return(result)
}

Lambda_new_generator = function(x, Theta_list, Z, eta_ij, Psi_ij)
{
  inverse_term = 0
  not_inv_term = 0
  
  for (i in 1:ncol(Z))
  {
    Theta_i = Theta_list[[i]]
    Sigma_inv = solve(Theta_i$Sigma)
    xi_i = matrix(Theta_i$xi, nrow = length(Theta_i$xi))
    
    for (j in 1:nrow(Z))
    {
      eta_ij_vec = t(eta_ij[[i]][j, , drop = F])
      x_ij = t(x[j, , drop = F])
      
      inverse_term = inverse_term + ( (Z[j, i]) * Psi_ij[[i]][[j]] * Sigma_inv)
      not_inv_term = not_inv_term + ( (Z[j, i]) * ( eta_ij_vec * (Sigma_inv %*% (x_ij - xi_i)) ) )
    }
  }
  
  Lambda_new = solve(inverse_term) %*% not_inv_term
  Lambda_new = diag(as.vector(Lambda_new), nrow = length(Lambda_new), ncol = length(Lambda_new))
  
  return(Lambda_new)
}

Lambda_generator_v2 = function(x, Theta_list)
{
  num_support = length(Theta_list)
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
  
  inv_not_inv_list =
    lapply(X = Theta_list,
           FUN = 
      function(Theta_i) 
      {
        n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
        p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
        
        Sigma_inv = solve(Theta_i$Sigma)
        xi_i_mat = matrix(Theta_i$xi, byrow = T, 
                          nrow = n, ncol = p)
        
        eta_i = Theta_i$eta_i
        Psi_i = Theta_i$Psi_i
        Z_i = Theta_i$Z_i
        
        inverse_term =
          apply(Psi_i, MARGIN = 3, FUN = function(Psi){Psi*Sigma_inv}, simplify = T) %>%
          apply(MARGIN = 1, FUN = function(Psi){sum(Psi*Z_i)}) %>%
          matrix(nrow = p, ncol = p)
        
        not_inv_term = 
          (Z_i * (eta_i * ((x - xi_i_mat) %*% Sigma_inv))) %>%
          apply(2, sum)
        
        result_mat = rbind(inverse_term, not_inv_term)
        
        return(result_mat)
      })
    
  inv_not_inv_mat =
    inv_not_inv_list %>%
    unlist() %>%
    array(dim = c(p+1, p, num_support)) %>%
    apply(MARGIN = c(1, 2), FUN = sum)
  
  inverse_term = inv_not_inv_mat[1:p, 1:p]
  not_inv_term = inv_not_inv_mat[p+1, ]
  
  Lambda_new = solve(inverse_term) %*% not_inv_term
  Lambda_new = diag(as.vector(Lambda_new), nrow = length(Lambda_new), ncol = length(Lambda_new))
  
  return(Lambda_new)
}

Lambda_generator_v2_lin = function(x, Theta_list)
{
  num_support = length(Theta_list)
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
  
  inv_not_inv_list =
    lapply(X = Theta_list,
           FUN = 
             function(Theta_i) 
             {
               n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
               p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
               
               Sigma_inv = solve(Theta_i$Sigma)
               xi_i_mat = matrix(Theta_i$xi, byrow = T, 
                                 nrow = n, ncol = p)
               
               eta_i = Theta_i$eta_i
               Psi_i = Theta_i$Psi_i
               Z_i = Theta_i$Z_i
               
               inverse_term =
                 Sigma_inv *
                 apply(Psi_i, MARGIN = c(1, 2), FUN = function(Psi){sum(Z_i * Psi)})
            
               eta_x_xi_mat = cbind(eta_i, (x - xi_i_mat))
              
               not_inv_term_Z_eta_x_xi =
                 eta_x_xi_mat %>%
                 apply(1, FUN = function(eta_x_xi_vec)
                   {outer(X = eta_x_xi_vec[1:p], Y = eta_x_xi_vec[(p+1):(2*p)])}) %>%
                 array(dim = c(p, p, n)) %>%
                 apply(MARGIN = c(1, 2), FUN = function(eta_x_xi_outer){sum(Z_i * eta_x_xi_outer)} )
                 
               not_inv_term = t((Sigma_inv * not_inv_term_Z_eta_x_xi) %*% rep(1, p))
               
               result_mat = rbind(inverse_term, not_inv_term)
               
               return(result_mat)
             })
  
  inv_not_inv_mat =
    inv_not_inv_list %>%
    unlist() %>%
    array(dim = c(p+1, p, num_support)) %>%
    apply(MARGIN = c(1, 2), FUN = sum)
  
  inverse_term = inv_not_inv_mat[1:p, 1:p]
  not_inv_term = inv_not_inv_mat[p+1, ]
  
  Lambda_new = solve(inverse_term) %*% not_inv_term
  Lambda_new = diag(as.vector(Lambda_new), nrow = length(Lambda_new), ncol = length(Lambda_new))
  
  return(Lambda_new)
}

Sigma_new_generator = function(x, Theta_list_new_Lambda, Q_support, Z, eta_ij, Psi_ij)
{
  sigma_new_numerator = 0
  
  for (i in 1:length(Theta_list_new_Lambda))
  {
    Theta_i = Theta_list_new_Lambda[[i]]
    alpha_i = Q_support$alpha[i]
    xi_i = Theta_i$xi
    Lambda_i = Theta_i$Lambda
    eta_i = eta_ij[[i]]
    
    for (j in 1:nrow(eta_i))
    {
      cross_term = (x[j,] - xi_i - (Lambda_i %*% eta_i[j,])) %>% as.numeric()
      first_term = (Z[j,i]/alpha_i) * (cross_term %*% t(cross_term))
      
      inside_term = Psi_ij[[i]][[j]] - (eta_i[j,] %*% t(eta_i[j,]))
      second_term = (Z[j,i]/alpha_i) * t(Lambda_i %*% inside_term %*% t(Lambda_i))
      
      sigma_new_numerator = sigma_new_numerator + (first_term + second_term)
    }
  }
  
  sigma_new = sigma_new_numerator/sum(Z)
  
  return(sigma_new)
}

Sigma_generator_v2 = function(x, Theta_list, Z)
{
  num_support = length(Theta_list)
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
  
  num_list = 
    lapply(X = Theta_list,
           FUN = 
    function(Theta_i)
    { 
      n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
      p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
      
      xi_i_mat = matrix(Theta_i$xi, byrow = T, 
                        nrow = n, ncol = p)
      Lambda_i = Theta_i$Lambda
      
      Z_i_div_alpha = Theta_i$Z_i / Theta_i$alpha_i
      Z_i_div_alpha[Z_i_div_alpha < 1e-32] = 0
      
      eta_i = Theta_i$eta_i 
      Psi_i = Theta_i$Psi_i
      
      cross_term = sqrt(Z_i_div_alpha) * (x - xi_i_mat - (eta_i %*% Lambda_i))
      first_term = t(cross_term) %*% cross_term
      
      eta_i_outer =
        apply(eta_i, 1, FUN = function(eta){eta %*% t(eta)}, simplify = T) %>%
        array(dim = c(p, p, n))
      
      second_term =
        (Psi_i - eta_i_outer) %>%
        apply(MARGIN = 3, FUN = function(inside_term){Lambda_i %*% inside_term %*% t(Lambda_i)}) %>%
        apply(MARGIN = 1, FUN = function(later_term){Z_i_div_alpha * later_term}) %>%
        apply(MARGIN = 2, FUN = sum) %>%
        matrix(nrow = p, ncol = p)
      
      result = first_term + second_term
      return(result)
    })
  
  Sigma_new_numerator =
    num_list %>%
    unlist() %>%
    array(dim = c(p, p, num_support)) %>%
    apply(MARGIN = c(1, 2), FUN = sum)
  
  Sigma_new = Sigma_new_numerator/sum(Z)
  
  return(Sigma_new)
}


Q_function = function(x, Theta_list)
{
  result_vec = 
    lapply(X = Theta_list,
           FUN = function(Theta_i)
           {
             n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
             p = ifelse(is.null(dim(x)), yes = 1, no = ncol(x))
             
             Z_i = Theta_i$Z_i
             eta_i = Theta_i$eta_i
             Psi_i = Theta_i$Psi_i
             
             w_i = Theta_i$prob
             xi_i_mat = 
               Theta_i$xi %>% 
               matrix(nrow = nrow(x),
                      ncol = ncol(x), byrow = T)
             Lambda_i = Theta_i$Lambda
             Sigma_inv = solve(Theta_i$Sigma)
             
             term_w = log(w_i)
             term_Sigma_inv = log(det(Sigma_inv))/2
             
             x_c = x - xi_i_mat - (eta_i %*% t(Lambda_i))
             
             term_cross_product = (((x_c %*% Sigma_inv)*x_c) %>% apply(1, sum))/(-2)
             
             eta_i_outer =
               apply(eta_i, 1, FUN = function(eta){eta %*% t(eta)}, simplify = T) %>%
               array(dim = c(p, p, n))
             inside_term = Psi_i - eta_i_outer
             
             term_trace = 
               (inside_term %>%
                  apply(MARGIN = 3, 
                        FUN = function(inside_i){ tr(Sigma_inv %*% Lambda_i %*% inside_i %*% t(Lambda_i)) }))/(-2)
             
             term_whole = Z_i * (term_w + term_Sigma_inv + term_cross_product + term_trace)
             result = sum(term_whole)
             
             return(result)
           }) %>% unlist()
  
  result = sum(result_vec)
  
  return(result)
}

#6. EM function in CNM_SSMMSN function
EM_part = function(x, Q_support, max_iter_EM, tol_EM, em_log_lik_check = F)
{
  ##############
  #   MLE(EM)  # 
  ##############
  
  iter_inner = 0
  log_lik_diff_inner = tol_EM + 1
  log_lik_set = density_mixture_v2(x, Q_support = Q_support) %>% log() %>% sum()
  
  if (em_log_lik_check == T)
  {
    result_prob_em = NULL
    result_xi_em = NULL
    result_Lambda_em = NULL
    result_Sigma_em = NULL
  }
  
  while ((iter_inner < max_iter_EM) & (log_lik_diff_inner > tol_EM)) 
  {
    iter_inner = iter_inner + 1
    
    ##########################
    # MLE(EM) - EM Starts #
    ##########################
    
    ##########################
    # MLE(EM) - 1. E-Step  #
    ##########################
    num_support = length(Q_support$alpha)
    
    Theta_list = Theta_list_generator(Q_support = Q_support)
    
    ##########################
    # MLE(EM) - 1-1. E(Z_ij) #
    ##########################
    ##########################
    # MLE(EM) - 1-2. E(tau) = eta and E(tau %*% t(tau)) = Psi #
    ##########################
    
    Z = Z_generator(x = x, Theta_list = Theta_list)
    mu_tau_TN =  lapply(X = Theta_list, FUN = mu_tau_TN_generator, x = x)
    eta_ij = lapply(X = mu_tau_TN, FUN = eta_generator)
    Psi_ij = lapply(X = mu_tau_TN, FUN = Psi_generator) #Psi_ij[[1]][[20]] i = 1, j = 20
    
    ##########################
    # MLE(EM) - 2. M-Step  #
    ##########################
    
    ##########################
    # MLE(EM) - 2-1. Update prob, w  #
    ##########################
    
    prob_new = apply(Z, 2, FUN = sum)/sum(Z)
    
    Q_support_new_prob = Q_support
    Q_support_new_prob$prob = prob_new
    
    ### Check whether Updating "Prob" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_prob = density_mixture_v2(x, Q_support = Q_support_new_prob) %>% log() %>% sum()
      result_prob_em = c(result_prob_em, Q_result_prob)
    }
    
    ##########################
    # MLE(EM) - 2-2. Update location parameter, xi  #
    ##########################
    Theta_list_new_prob = 
      Theta_list_generator(Q_support = Q_support_new_prob, 
                           Z = Z, eta_ij = eta_ij, Psi_ij = Psi_ij)
    
    xi_new = xi_generator_v2(Theta_list = Theta_list_new_prob, x = x)
    
    Q_support_new_xi = Q_support_new_prob
    Q_support_new_xi$xi = xi_new
    
    ### Check whether Updating "Xi" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_xi = density_mixture_v2(x, Q_support = Q_support_new_xi) %>% log() %>% sum()
      result_xi_em = c(result_xi_em, Q_result_xi)
    }
    
    ##########################
    # MLE(EM) - 2-3. Update skewness parameter, Lambda#
    ##########################
    Theta_list_new_xi = 
      Theta_list_generator(Q_support = Q_support_new_xi, 
                           Z = Z, eta_ij = eta_ij, Psi_ij = Psi_ij)
  
    Lambda_new =
      Lambda_generator_v2(x = x, Theta_list = Theta_list_new_xi)
    
    Q_support_new_Lambda = Q_support_new_xi
    Q_support_new_Lambda$Lambda = Lambda_new
    
    ### Check whether Updating "Lambda" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Lambda = density_mixture_v2(x, Q_support = Q_support_new_Lambda) %>% log() %>% sum()
      result_Lambda_em = c(result_Lambda_em, Q_result_Lambda)
    }
    
    ##########################
    # MLE(EM) - 2-4. Update scale parameter, Sigma#
    ##########################
    Theta_list_new_Lambda = 
      Theta_list_generator(Q_support = Q_support_new_Lambda, 
                           Z = Z, eta_ij = eta_ij, Psi_ij = Psi_ij)
    
    Sigma_new = Sigma_generator_v2(x = x,
                                   Theta_list = Theta_list_new_Lambda, Z = Z)
    
    Q_support_new_Sigma = Q_support_new_Lambda
    Q_support_new_Sigma$Sigma = (Sigma_new + t(Sigma_new))/2
    
    ### Check whether Updating "Sigma" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Sigma = density_mixture_v2(x, Q_support = Q_support_new_Sigma) %>% log() %>% sum()
      result_Sigma_em = c(result_Sigma_em, Q_result_Sigma)
    }
    
    ##########################
    # MLE(EM) - 2-5. Save all updated parameters
    ##########################
    Q_support = Q_support_new_Sigma
    
    ##########################
    # MLE(EM) - 2-6. Calculate to identify wheter stopping rules are satisified or not
    ##########################
    log_lik_new_em = density_mixture_v2(x, Q_support = Q_support) %>% log() %>% sum()
    
    log_lik_set = c(log_lik_set, log_lik_new_em)
    
    log_lik_diff_inner = log_lik_set[(iter_inner+1)] - log_lik_set[iter_inner] 
  }
  
  ##########################
  # MLE(EM) - 2-7 Save all updated parameters after EM-algorithm is done
  ##########################
  result_list = list()
  
  result_list$Q_support = Q_support
  
  result_list$iter_inner = iter_inner
  
  if(em_log_lik_check == T)
  {
    result_list$result_prob_em = result_prob_em
    result_list$result_xi_em = result_xi_em 
    result_list$result_Lambda_em = result_Lambda_em
    result_list$result_Sigma_em = result_Sigma_em
  }
  
  return(result_list)
}

central_difference = function(func, value, h = 1e-4, ...)
{
  if ((length(value) != 1)) 
  {
    print("value should be single numeric")
    return(NULL)
  }
  numerator = func(value + h, ...) - func(value - h, ...)
  denominator = 2*h
  
  result = numerator/denominator
  
  return(result)
}

where_potential_maximum_grid = function(logical_vector)
{
  target_vec = logical_vector
  len_vec = length(target_vec)
  
  A_vec = target_vec[-len_vec]
  B_vec = target_vec[-1]
  
  C_vec = (A_vec == T) & (B_vec == F)
  
  start_ind_2 = (1:(len_vec-1))[C_vec]
  end_ind_2 = start_ind_2 + 1
  
  if (target_vec[1] == F) {
    start_ind_1 = 1
    end_ind_1 = 2
  } else {
    start_ind_1 = NULL
    end_ind_1 = NULL
  }
  
  if (target_vec[len_vec] == T) {
    start_ind_3 = len_vec - 1
    end_ind_3 = len_vec
  } else {
    start_ind_3 = NULL
    end_ind_3 = NULL
  }
  
  start_ind = c(start_ind_1, start_ind_2, start_ind_3)
  end_ind = c(end_ind_1, end_ind_2, end_ind_3)
  
  mat_ind = cbind(start_ind, end_ind)
  colnames(mat_ind) = c("start", "end")
  
  return(mat_ind)
}

give_interval_with_potential = function(scale_interval, grid_length, func, 
                                        x, Q_support, 
                                        maximum_grid_length = 2*grid_length)
{
  grid_for_gradient = seq(from = min(scale_interval), 
                          to = max(scale_interval),
                          length.out = grid_length)
  
  derivatives_of_gradient_function =
    sapply(X = grid_for_gradient, 
           FUN = function(value){
             central_difference(func = func, 
                                value = value, Q_support = Q_support, x = x)
           })
  
  derivatives_positive = derivatives_of_gradient_function > 0
  
  grid_length_temp = grid_length
  scale_interval_temp = scale_interval
  
  inc_grid_length_by = ceiling(grid_length/5)
  
  while ( (rev(derivatives_positive)[1] == T) & (length(derivatives_positive) < maximum_grid_length) )
  {
    grid_gap = grid_for_gradient[2] - grid_for_gradient[1]
    
    scale_interval_temp[which.max(scale_interval_temp)] = max(scale_interval_temp) + inc_grid_length_by*grid_gap 
    
    grid_length_old = grid_length_temp
    grid_length_temp = grid_length_temp + inc_grid_length_by
    
    grid_for_gradient = seq(from = min(scale_interval_temp), 
                            to = max(scale_interval_temp),
                            length.out = grid_length_temp) 
    
    additional_grid =
      grid_for_gradient[(grid_length_old+1):(grid_length_old + inc_grid_length_by)]
    
    additional_derivative =
      sapply(X = additional_grid, 
             FUN = function(value){
               central_difference(func = func, 
                                  value = value, Q_support = Q_support, x = x)
             })
    
    derivatives_positive = c(derivatives_positive, (additional_derivative > 0))
  }
  
  grid_with_potential_mat = where_potential_maximum_grid(logical_vector = derivatives_positive)
  
  result_list = list()
  result_list$grid_for_gradient = grid_for_gradient
  result_list$grid_with_potential_mat = grid_with_potential_mat
  
  return(result_list)
}

#5. CNM function in CNM_SSMMSN function 

CNM_part = function(x, Q_support, scale_interval, grid_length, 
                    gamma_tol = length(x)/length(Q_support$xi) * 1e-14)
{
  scale_interval_free = 
    c(min(scale_interval), 
      max(scale_interval, Q_support$alpha))
  
  if (length(Q_support$alpha) == 1) #Do VEM when only one support exists
  {
    interval_with_potential_list =
      give_interval_with_potential(scale_interval = scale_interval_free, grid_length = grid_length,
                                   func = Directional, 
                                   x = x, Q_support = Q_support)
    
    grid_with_potential_mat = interval_with_potential_list$grid_with_potential_mat
    grid_for_gradient = interval_with_potential_list$grid_for_gradient
    
    gradient_results = 
      apply(X = grid_with_potential_mat, MARGIN = 1,
            FUN = function(index){
              optimize(Directional, 
                       interval = c(grid_for_gradient[index[1]], grid_for_gradient[index[2]]), 
                       maximum = TRUE, x = x, Q_support = Q_support)
            })
    
    gradient_results_mat = matrix(unlist(gradient_results), ncol = 2, byrow = T)
    colnames(gradient_results_mat) = c("maximum", "objective")
    
    idx_maximizer_among_grids = which.max(gradient_results_mat[,"objective"])
    
    Directional_Derivative = gradient_results_mat[idx_maximizer_among_grids, "objective"]
    maximizer_DD = gradient_results_mat[idx_maximizer_among_grids, "maximum"]
    
    DD_with_Q_support =
      sapply(X = Q_support$alpha, FUN = Directional, 
             x = x, Q_support = Q_support)
    
    minimizer_DD_idx = which.min(DD_with_Q_support)
    
    optimized_L = optimize(L_function_vem, interval = c(0,1), maximum = T,
                           x = x, Q_support = Q_support,
                           maximizer_parameter = maximizer_DD, minimizer_idx = minimizer_DD_idx)
    
    maximizer_L_alpha = optimized_L$maximum
    
    ## Update
    minizer_prob_in_support = Q_support$prob[minimizer_DD_idx]
    
    Q_support$prob[minimizer_DD_idx] = minizer_prob_in_support * (1 - maximizer_L_alpha)
    
    Q_support$alpha = c(Q_support$alpha, maximizer_DD)
    Q_support$prob = c(Q_support$prob, 
                       minizer_prob_in_support*maximizer_L_alpha)
  } else { #Do CNM
    
    interval_with_potential_list =
      give_interval_with_potential(scale_interval = scale_interval_free, grid_length = grid_length,
                                   func = Directional_with_identifiability, 
                                   x = x, Q_support = Q_support)
    
    grid_with_potential_mat = interval_with_potential_list$grid_with_potential_mat
    grid_for_gradient = interval_with_potential_list$grid_for_gradient
    
    gradient_results = 
      apply(X = grid_with_potential_mat, MARGIN = 1,
            FUN = function(index){
              optimize(Directional_with_identifiability, 
                       interval = c(grid_for_gradient[index[1]], grid_for_gradient[index[2]]), 
                       maximum = TRUE, x = x, Q_support = Q_support)
            })
    
    gradient_results_mat = matrix(unlist(gradient_results), ncol = 2, byrow = T)
    colnames(gradient_results_mat) = c("maximum", "objective")
    
    idx_maximizer_among_grids = which.max(gradient_results_mat[,"objective"])
    
    Directional_Derivative_set = gradient_results_mat[, "objective"]
    maximizer_DD_set = gradient_results_mat[, "maximum"]
    
    Directional_Derivative = Directional_Derivative_set[idx_maximizer_among_grids]
    
    Q_support_plus = Q_support
    
    Q_support_plus$alpha = c(Q_support$alpha, maximizer_DD_set)
    
    Q_support_plus$prob = c(Q_support$prob, 
                            rep(0, times = (length(Q_support_plus$alpha) - length(Q_support$alpha))))
    
    S_denominator_vector = density_mixture_v2(x, Q_support = Q_support_plus)
    S_denominator = matrix(S_denominator_vector, nrow = length(x)/length(Q_support$xi), ncol = length(Q_support_plus$alpha))
    
    S_numerator = matrix(data = 0, nrow = nrow(S_denominator), ncol = ncol(S_denominator))
    
    for (i in 1:length(Q_support_plus$alpha))
    {
      S_numerator_col_i =
        dmsn_sahu_v3_biv(x = x, 
                         xi = Q_support_plus$xi, Lambda = Q_support_plus$Lambda, 
                         Sigma = Q_support_plus$alpha[i]*Q_support_plus$Sigma)
      
      S_numerator[,i] = S_numerator_col_i
    }
    
    S_mat = S_numerator/S_denominator
    
    x_4nnls = rbind(S_mat*sqrt(gamma_tol), rep(1, ncol(S_mat)))
    y_4nnls = c(rep(2,nrow(S_mat))*sqrt(gamma_tol), 1)
    
    nnls_result = nnls(x_4nnls, y_4nnls)
    Q_support_plus$prob = nnls_result$x
    
    idx_not_zero = nnls_result$x != 0
    idx_not_zero[1] = TRUE
    
    Q_support$prob = Q_support_plus$prob[idx_not_zero]
    Q_support$prob = Q_support$prob / sum(Q_support$prob)
    Q_support$alpha = Q_support_plus$alpha[idx_not_zero]
  }
  
  result_list = list()
  result_list$Directional_Derivative = Directional_Derivative
  result_list$alpha = Q_support$alpha
  result_list$prob = Q_support$prob
  
  return(result_list)
}

CNM_SSMMSN = function(x, Q_support_initial,
                      tol_outer = 1, tol_Directional_Derivative = 1e-6, 
                      tol_EM = 0.1,
                      max_iter_outer = 10, max_iter_EM = 1,
                      scale_interval = c(1, 1000), grid_length = 100,
                      print_each_iteration = T)
{
  Q_support = Q_support_initial
  
  log_lik_diff_outer = tol_outer + 1
  Directional_Derivative = tol_Directional_Derivative + 1
  iter = 0

  LL_set = density_mixture_v2(x, Q_support = Q_support_initial) %>% log() %>% sum()
  DD_set = NULL
  
  while((iter < max_iter_outer) & (log_lik_diff_outer > tol_outer)) # 
  {
    start_time = Sys.time()
    iter = iter + 1
    
    Q_support_old = Q_support
    
    ##############
    # NPMLE(CNM) #
    ##############
    
    CNM_result_list = 
      CNM_part(x = x, Q_support = Q_support, 
               scale_interval = scale_interval, grid_length = grid_length)
    
    Directional_Derivative = CNM_result_list$Directional_Derivative
    
    if (Directional_Derivative >= tol_Directional_Derivative)
    {
      Q_support$alpha = CNM_result_list$alpha
      Q_support$prob = CNM_result_list$prob
    } else if (print_each_iteration == T) {
      cat("Maximum Gradient: ", Directional_Derivative,  "under tolerance: ", tol_Directional_Derivative, "\n")
      cat("Therefore, no update with CNM \n")
    }
    
    L_inter = density_mixture_v2(x, Q_support = Q_support) %>% log() %>% sum()
    
    ##############
    #   MLE(EM)  # Find Fixed Parameters with EM algorithm
    ##############
    
    EM_result_list = 
      EM_part(x = x, 
              Q_support = Q_support,
              tol_EM = tol_EM, max_iter_EM = max_iter_EM)    
    
    Q_support_new = EM_result_list$Q_support
    
    Q_support = Q_support_new
    
    ##############
    #   Calculation for stopping criterion & Printing result
    ##############
    
    Log_Likelihood_old = density_mixture_v2(x, Q_support = Q_support_old) %>% log() %>% sum()
    Log_Likelihood = density_mixture_v2(x, Q_support = Q_support) %>% log() %>% sum()
    
    spent_time = Sys.time() - start_time
    
    LL_set = c(LL_set, Log_Likelihood)
    DD_set = c(DD_set, Directional_Derivative)
    
    log_lik_diff_outer = LL_set[iter + 1] - LL_set[iter]
    
    if (print_each_iteration == T)
    {
      cat("########################################################", "\n")
      cat("iteration spent for EM = ", EM_result_list$iter_inner, "\n")
      cat("iteration = " , iter , 
          "D =", Directional_Derivative,  
          "L = ", Log_Likelihood, 
          "L_inter = ", L_inter,
          "L_old = ", Log_Likelihood_old, 
          "Support =", length(Q_support$alpha), "\n")
      cat(Q_support$alpha, "\n")
      cat("Spent time: ") 
      print(spent_time) 
    }
  }
  
  ##############
  # Save the resulting parameters and return the result
  ##############
  
  result_list = list()
  result_list$Q_support = Q_support
  result_list$Log_likelihood = LL_set
  result_list$Directonal_derivative = DD_set
  
  return(result_list)
}

library(cubature)

ISE_adaptintegral = function(f_true, f_est, x_set)
{
  f_temp = function(x) {(f_true(x) - f_est(x))^2}
  lower_temp = x_set %>% apply(MARGIN = 2, FUN = min)
  upper_temp = x_set %>% apply(MARGIN = 2, FUN = max)
  
  result_integrate = adaptIntegrate(f = f_temp, lowerLimit = lower_temp, upperLimit = upper_temp)
  
  return(result_integrate)
}


ISE_MC = function(f_true, f_est, X_sample)
{
  fn = f_est(x = X_sample) 
  f0 = f_true(x = X_sample)
  
  result = mean(((fn - f0)^2)/f0, na.rm = T)
  
  return(result)
}

MSN_initial_generator = function(x)
{
  x_bar = x %>% apply(MARGIN = 2, mean)
  x_bar_mat = x_bar %>% matrix(nrow = nrow(x), ncol = ncol(x), byrow = T)
  
  x_c = x - x_bar_mat
  
  sample_skewness = apply((x_c^3), 2, mean)/(apply((x_c^2), 2, mean)^(3/2))
  sign_sample_skewness = sample_skewness/abs(sample_skewness)
  
  Cov_mat = (1/nrow(x))*(t(x_c) %*% (x_c))
  I_mat = diag(nrow = nrow(Cov_mat))

  Cov_sqrt_diag_inv = diag(1/sqrt(diag(Cov_mat)), nrow = nrow(Cov_mat))
  
  Corr_initial = Cov_sqrt_diag_inv %*% Cov_mat %*% Cov_sqrt_diag_inv
  
  min_a = min(abs(Corr_initial))
  
  optimized_a =
    optimize(f = function(a)
                  {
                    Sigma_initial = Cov_mat - (1 - a)*(Cov_mat*I_mat)
                    
                    Lambda_initial = sign_sample_skewness * sqrt( (pi*(1-a)/(pi - 2)) * (Cov_mat*I_mat))
                    
                    xi_initial = x_bar - (sqrt(2/pi) * diag(Lambda_initial))
                    
                    Q_support_temp = list()
                    Q_support_temp$Sigma = Sigma_initial
                    Q_support_temp$Lambda = Lambda_initial
                    Q_support_temp$xi = xi_initial
                    Q_support_temp$alpha = 1
                    Q_support_temp$prob = 1
                    
                    log_lik_temp = density_mixture_v2(x = x, Q_support = Q_support_temp) %>% log() %>% sum()
                    
                    return(log_lik_temp)
                  }, maximum = T, interval = c(min_a, 1))
  
  Sigma_initial = Cov_mat - (1 - optimized_a$maximum)*(Cov_mat*I_mat)
  
  Lambda_initial = sign_sample_skewness * sqrt( (pi*(1 - optimized_a$maximum)/(pi - 2)) * (Cov_mat*I_mat))
  
  xi_initial = x_bar - (sqrt(2/pi) * diag(Lambda_initial))
  
  result_list = list()
  result_list$Sigma_init = Sigma_initial
  result_list$Lambda_init = Lambda_initial
  result_list$xi_init = xi_initial
  
  return(result_list)
}

MSN_initial_zero_Lambda_generator = function(x)
{
  x_bar = x %>% apply(MARGIN = 2, mean)
  x_bar_mat = x_bar %>% matrix(nrow = nrow(x), ncol = ncol(x), byrow = T)
  
  x_c = x - x_bar_mat
  
  sample_skewness = apply((x_c^3), 2, mean)/(apply((x_c^2), 2, mean)^(3/2))
  sign_sample_skewness = sample_skewness/abs(sample_skewness)
  
  Cov_mat = (1/nrow(x))*(t(x_c) %*% (x_c))
  I_mat = diag(nrow = nrow(Cov_mat))
  
  Cov_sqrt_diag_inv = diag(1/sqrt(diag(Cov_mat)), nrow = nrow(Cov_mat))
  
  Corr_initial = Cov_sqrt_diag_inv %*% Cov_mat %*% Cov_sqrt_diag_inv
  
  # min_a = min(abs(Corr_initial))
  # a = max(min_a, 0.95)
  
  # Sigma_initial = Cov_mat - (1 - a)*(Cov_mat*I_mat)
  Sigma_initial = Cov_mat
  
  # Lambda_initial = sign_sample_skewness * sqrt( (pi*(1 - optimized_a$maximum)/(pi - 2)) * (Cov_mat*I_mat))
  # Lambda_initial = sign_sample_skewness * sqrt( (pi*(1 - a)/(pi - 2)) * (Cov_mat*I_mat)) 
  Lambda_initial = 0 * (Cov_mat*I_mat) # If you set optimized_a$maximum = 1, then you will get this value
  
  xi_initial = x_bar - (sqrt(2/pi) * diag(Lambda_initial))
  
  result_list = list()
  result_list$Sigma_init = Sigma_initial
  result_list$Lambda_init = Lambda_initial
  result_list$xi_init = xi_initial
  
  return(result_list)
}

# Mixture of Mixture 220710

FMSN_initial_generator = function(x, L, initial_clustering_attempt = 4, 
                                  how_much_initial = 5,
                                  label_count_minimum = 5)
{
  clustering_attempt_per_method = round(initial_clustering_attempt/2)
  
  label_est_kmeans_list = list()
  
  for (kmeans_attempt_i in 1:clustering_attempt_per_method)
  {
    kmeans_x_4_initial = kmeans(x = x, centers = L, nstart = 1)
    
    label_est_kmeans_list[[kmeans_attempt_i]] = kmeans_x_4_initial$cluster
  }
  
  label_est_gmm_list = list()
  
  suppressPackageStartupMessages(expr = library(mclust))
  for (gmm_attempt_i in 1:clustering_attempt_per_method)
  {
    gmm_x_4_initial = mclust::Mclust(data = x, G = L, verbose = F)
    
    label_est_gmm_list[[gmm_attempt_i]] = gmm_x_4_initial$classification
  }
  detach("package:mclust", unload=TRUE)
  
  label_est_list = append(label_est_kmeans_list, label_est_gmm_list)
  
  unvalid_label_idx = 
    label_est_list %>%
    lapply(FUN = function(label_est){
      min_label = label_est %>% table() %>% min()
      is_label_unvalid = (min_label < label_count_minimum)
    }) %>% unlist()
  
  label_est_list = label_est_list[!unvalid_label_idx]
  label_est_list = label_est_list[!duplicated(label_est_list)]
  
  result_per_attempt = 
    lapply(X = label_est_list,
           FUN = function(label_est)
           {
             Q_support_initial_list = list()
             
             for(l in 1:L)
             {
               x_l_temp = x[label_est == l, , drop = F]
               
               x_l_initial = MSN_initial_generator(x = x_l_temp)
               
               Q_support_l_initial = list()
               Q_support_l_initial$alpha = c(1)
               Q_support_l_initial$xi = x_l_initial$xi_init
               Q_support_l_initial$Sigma = x_l_initial$Sigma_init
               Q_support_l_initial$Lambda = x_l_initial$Lambda_init
               Q_support_l_initial$prob = (label_est == l) %>% mean()
               
               Q_support_initial_list[[l]] = Q_support_l_initial
             }
             
             log_lik_temp = density_mixture_Q_support_list(x = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
             
             result_list = list()
             result_list$Q_support_initial_list = Q_support_initial_list
             result_list$log_lik = log_lik_temp
             
             return(result_list)
           })
  
  log_lik_set = rep(NA, times = length(result_per_attempt))
  for (attempt_i_valid in 1:length(log_lik_set))
  {
    log_lik_set[attempt_i_valid] = 
      result_per_attempt[[attempt_i_valid]]$log_lik
  }
  
  unique_log_lik_idx = !duplicated(log_lik_set)
  unique_log_lik_order = unique_log_lik_idx*order(log_lik_set, decreasing = T)
  
  unique_log_lik_order_under_valid = (unique_log_lik_order[unique_log_lik_idx] %>% order()) <= how_much_initial
  
  selected_idx = (unique_log_lik_order[unique_log_lik_idx])[unique_log_lik_order_under_valid]
  
  selected_list_number = (1:length(log_lik_set))[unique_log_lik_order %in% selected_idx]
  
  Q_support_initial_list_list = list()
  
  for (i in 1:length(selected_list_number))
  {
    Q_support_initial_list_list[[i]] = 
      result_per_attempt[[ (selected_list_number[i]) ]]$Q_support_initial_list
  }
  
  return(Q_support_initial_list_list)
}

FMSN_initial_generator_true = function(x, label_true)
{
   Q_support_initial_list = list()
   
   unique_label = unique(label_true)
   
   for(l in 1:length(unique_label))
   {
     x_l_temp = x[label_true == unique_label[l], , drop = F]
     
     x_l_initial = MSN_initial_generator(x = x_l_temp)
     
     Q_support_l_initial = list()
     Q_support_l_initial$alpha = c(1)
     Q_support_l_initial$xi = x_l_initial$xi_init
     Q_support_l_initial$Sigma = x_l_initial$Sigma_init
     Q_support_l_initial$Lambda = x_l_initial$Lambda_init
     Q_support_l_initial$prob = (label_true == unique_label[l]) %>% mean()
     
     Q_support_initial_list[[l]] = Q_support_l_initial
   }
   
  log_lik_temp = density_mixture_Q_support_list(x = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
  
  return(Q_support_initial_list)
}

fmmst_initial_generator = function(Q_support_initial_list)
{
  L = length(Q_support_initial_list)
  
  fmmst_initial = list()
  fmmst_initial$mu = list()
  fmmst_initial$sigma = list()
  fmmst_initial$delta = list()
  fmmst_initial$pro = rep(NA, times = L)
  
  for (l in 1:L)
  {
    fmmst_initial$mu[[l]] = Q_support_initial_list[[l]]$xi %>% as.matrix()
    fmmst_initial$sigma[[l]] = Q_support_initial_list[[l]]$Sigma 
    fmmst_initial$delta[[l]] = diag(Q_support_initial_list[[l]]$Lambda) %>% as.matrix() 
    fmmst_initial$pro[l] =  Q_support_initial_list[[l]]$prob
  }
  
  fmmst_initial$dof = rep(4, times = L)
  
  return(fmmst_initial)
}

FMVN_initial_generator = function(Q_support_initial_list)
{
  prob_vec = 
    lapply(X = Q_support_initial_list,
           FUN = function(Q_support){
             Q_support$prob
           }) %>% unlist()
  mu_list =
    lapply(X = Q_support_initial_list,
           FUN = function(Q_support){
             Q_support$xi
           })
  Sigma_list = 
    lapply(X = Q_support_initial_list,
           FUN = function(Q_support){
             Q_support$Sigma + (1 - 2/pi)*(Q_support$Lambda^2)
           })
  
  result_list = list()
  result_list$prob_vec = prob_vec
  result_list$mu_list = mu_list
  result_list$Sigma_list = Sigma_list
  
  return(result_list)
}


FMSN_initial_5clustering = function(x, L, 
                                    how_much_initial = 5,
                                    how_much_per_clustering = 4, 
                                    label_count_minimum = 5)
{
  
  label_est_kmeans_list = list()
  for (kmeans_attempt_i in 1:how_much_per_clustering)
  {
    kmeans_x_4_initial = kmeans(x = x, centers = L, nstart = 20)
    label_est_kmeans_list[[kmeans_attempt_i]] = kmeans_x_4_initial$cluster
  }
  
  suppressPackageStartupMessages(expr = library(mclust))
  label_est_mhc_list = list()
  for (mhc_attempt_i in 1:how_much_per_clustering)
  {
    mhc_x_4_initial = hc(x) %>% hclass(G = L) %>% as.vector()
    label_est_mhc_list[[mhc_attempt_i]] = mhc_x_4_initial
  }
  
  label_est_gmm_list = list()
  for (gmm_attempt_i in 1:how_much_per_clustering)
  {
    gmm_x_4_initial = mclust::Mclust(data = x, G = L, verbose = F)
    
    label_est_gmm_list[[gmm_attempt_i]] = gmm_x_4_initial$classification
  }
  detach("package:mclust", unload=TRUE)
  
  label_est_hc_list = list()
  for (hc_attempt_i in 1:how_much_per_clustering)
  {
    hc_dist = dist(x = x)
    hc_x_4_initial = hclust(d = hc_dist, method = "ward.D2")
    label_est_hc_list[[hc_attempt_i]] = hc_x_4_initial %>% cutree(k = L)
  }
  
  label_est_specc_list = list()
  for (specc_attempt_i in 1:how_much_per_clustering)
  {
    if (L == 1)
    {
      label_est_specc_list[[specc_attempt_i]] = rep(1, nrow(x))
    } else {
      tryCatch(expr = {
        specc_x_4_initial = kernlab::specc(x = x, centers = L)
        label_est_specc_list[[specc_attempt_i]] = specc_x_4_initial@.Data
      }, error = function(e){print(e)})
    }
  }
  
  label_est_list = list()
  label_est_list = 
    label_est_list %>% 
    append(label_est_kmeans_list) %>%
    append(label_est_mhc_list) %>%
    append(label_est_gmm_list) %>%
    append(label_est_hc_list) %>%
    append(label_est_specc_list)
  
  null_label_idx = 
    label_est_list %>%
    lapply(FUN = function(label_est){
      is.null(label_est) %>% return()
    }) %>% unlist()
  
  label_est_list = label_est_list[!null_label_idx]
  
  unvalid_label_idx = 
    label_est_list %>%
    lapply(FUN = function(label_est){
      min_label = label_est %>% table() %>% min()
      is_label_unvalid = (min_label < label_count_minimum)
    }) %>% unlist()
  
  label_est_list = label_est_list[!unvalid_label_idx]
  label_est_list = label_est_list[!duplicated(label_est_list)]
  
  result_per_attempt = 
    lapply(X = label_est_list,
           FUN = function(label_est)
           {
             Q_support_initial_list = list()
             
             unique_label = unique(label_est)
             
             for(l in 1:length(unique_label))
             {
               x_l_temp = x[label_est == (unique_label[l]), , drop = F]
               
               x_l_initial = MSN_initial_generator(x = x_l_temp)
               
               Q_support_l_initial = list()
               Q_support_l_initial$alpha = c(1)
               Q_support_l_initial$xi = x_l_initial$xi_init
               Q_support_l_initial$Sigma = x_l_initial$Sigma_init
               Q_support_l_initial$Lambda = x_l_initial$Lambda_init
               Q_support_l_initial$prob = (label_est == (unique_label[l])) %>% mean()
               
               Q_support_initial_list[[l]] = Q_support_l_initial
             }
             
             log_lik_temp = density_mixture_Q_support_list(x  = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
             
             result_list = list()
             result_list$Q_support_initial_list = Q_support_initial_list
             result_list$log_lik = log_lik_temp
             
             return(result_list)
           })
  
  log_lik_set = rep(NA, times = length(result_per_attempt))
  for (attempt_i_valid in 1:length(log_lik_set))
  {
    log_lik_set[attempt_i_valid] = 
      result_per_attempt[[attempt_i_valid]]$log_lik
  }
  
  unique_log_lik_idx = !duplicated(log_lik_set)
  unique_log_lik_order = unique_log_lik_idx*order(log_lik_set, decreasing = T)
  
  unique_log_lik_order_under_valid = (unique_log_lik_order[unique_log_lik_idx] %>% order()) <= how_much_initial
  
  selected_idx = (unique_log_lik_order[unique_log_lik_idx])[unique_log_lik_order_under_valid]
  
  selected_list_number = (1:length(log_lik_set))[unique_log_lik_order %in% selected_idx]
  
  Q_support_initial_list_list = list()
  
  for (i in 1:length(selected_list_number))
  {
    Q_support_initial_list_list[[i]] = 
      result_per_attempt[[ (selected_list_number[i]) ]]$Q_support_initial_list
  }
  
  return(Q_support_initial_list_list)
}


FMSN_initial_with_zero_Lambda_5clustering = 
  function(x, L, 
           how_much_initial = 10,
           how_much_per_clustering = 1, 
           label_count_minimum = 5)
{
    
    label_est_kmeans_list = list()
    for (kmeans_attempt_i in 1:how_much_per_clustering)
    {
      kmeans_x_4_initial = kmeans(x = x, centers = L, nstart = 20)
      label_est_kmeans_list[[kmeans_attempt_i]] = kmeans_x_4_initial$cluster
    }
    
    suppressPackageStartupMessages(expr = library(mclust))
    label_est_mhc_list = list()
    for (mhc_attempt_i in 1:how_much_per_clustering)
    {
      mhc_x_4_initial = hc(x) %>% hclass(G = L) %>% as.vector()
      label_est_mhc_list[[mhc_attempt_i]] = mhc_x_4_initial
    }
    
    label_est_gmm_list = list()
    for (gmm_attempt_i in 1:how_much_per_clustering)
    {
      gmm_x_4_initial = mclust::Mclust(data = x, G = L, verbose = F)
      
      label_est_gmm_list[[gmm_attempt_i]] = gmm_x_4_initial$classification
    }
    detach("package:mclust", unload=TRUE)
    
    label_est_hc_list = list()
    for (hc_attempt_i in 1:how_much_per_clustering)
    {
      hc_dist = dist(x = x)
      hc_x_4_initial = hclust(d = hc_dist, method = "ward.D2")
      label_est_hc_list[[hc_attempt_i]] = hc_x_4_initial %>% cutree(k = L)
    }
    
    label_est_specc_list = list()
    for (specc_attempt_i in 1:how_much_per_clustering)
    {
      if (L == 1)
      {
        label_est_specc_list[[specc_attempt_i]] = rep(1, nrow(x))
      } else {
        tryCatch(expr = {
          specc_x_4_initial = kernlab::specc(x = x, centers = L)
          label_est_specc_list[[specc_attempt_i]] = specc_x_4_initial@.Data
        }, error = function(e){print(e)})
      }
    }
    
    label_est_list = list()
    label_est_list = 
      label_est_list %>% 
      append(label_est_kmeans_list) %>%
      append(label_est_mhc_list) %>%
      append(label_est_gmm_list) %>%
      append(label_est_hc_list) %>%
      append(label_est_specc_list)
    
    null_label_idx = 
      label_est_list %>%
      lapply(FUN = function(label_est){
        is.null(label_est) %>% return()
      }) %>% unlist()
    
    label_est_list = label_est_list[!null_label_idx]
    
    unvalid_label_idx = 
      label_est_list %>%
      lapply(FUN = function(label_est){
        min_label = label_est %>% table() %>% min()
        is_label_unvalid = (min_label < label_count_minimum)
      }) %>% unlist()
    
    label_est_list = label_est_list[!unvalid_label_idx]
    label_est_list = label_est_list[!duplicated(label_est_list)]
    
    result_per_attempt = 
      lapply(X = label_est_list,
             FUN = function(label_est)
             {
               Q_support_initial_list = list()
               
               unique_label = unique(label_est)
               
               for(l in 1:length(unique_label))
               {
                 x_l_temp = x[label_est == (unique_label[l]), , drop = F]
                 
                 x_l_initial = MSN_initial_generator(x = x_l_temp)
                 
                 Q_support_l_initial = list()
                 Q_support_l_initial$alpha = c(1)
                 Q_support_l_initial$xi = x_l_initial$xi_init
                 Q_support_l_initial$Sigma = x_l_initial$Sigma_init
                 Q_support_l_initial$Lambda = x_l_initial$Lambda_init
                 Q_support_l_initial$prob = (label_est == (unique_label[l])) %>% mean()
                 
                 Q_support_initial_list[[l]] = Q_support_l_initial
               }
               
               log_lik_temp = density_mixture_Q_support_list(x = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
               
               result_list = list()
               result_list$Q_support_initial_list = Q_support_initial_list
               result_list$log_lik = log_lik_temp
               
               return(result_list)
             })
    
    result_with_zero_per_attempt = 
      lapply(X = label_est_list,
             FUN = function(label_est)
             {
               Q_support_initial_list = list()
               
               unique_label = unique(label_est)
               
               for(l in 1:length(unique_label))
               {
                 x_l_temp = x[label_est == (unique_label[l]), , drop = F]
                 
                 x_l_initial = MSN_initial_zero_Lambda_generator(x = x_l_temp)
                 
                 Q_support_l_initial = list()
                 Q_support_l_initial$alpha = c(1)
                 Q_support_l_initial$xi = x_l_initial$xi_init
                 Q_support_l_initial$Sigma = x_l_initial$Sigma_init
                 Q_support_l_initial$Lambda = x_l_initial$Lambda_init
                 Q_support_l_initial$prob = (label_est == (unique_label[l])) %>% mean()
                 
                 Q_support_initial_list[[l]] = Q_support_l_initial
               }
               
               log_lik_temp = density_mixture_Q_support_list(x = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
               
               result_list = list()
               result_list$Q_support_initial_list = Q_support_initial_list
               result_list$log_lik = log_lik_temp
               
               return(result_list)
             })
    
    result_per_attempt_whole = c(result_per_attempt, result_with_zero_per_attempt)
    
    log_lik_set = rep(NA, times = length(result_per_attempt_whole))
    for (attempt_i_valid in 1:length(log_lik_set))
    {
      log_lik_set[attempt_i_valid] = 
        result_per_attempt_whole[[attempt_i_valid]]$log_lik
    }
    
    unique_log_lik_idx = !duplicated(log_lik_set)
    unique_log_lik_order = unique_log_lik_idx*order(log_lik_set, decreasing = T)
    
    unique_log_lik_order_under_valid = (unique_log_lik_order[unique_log_lik_idx] %>% order()) <= how_much_initial
    
    selected_idx = (unique_log_lik_order[unique_log_lik_idx])[unique_log_lik_order_under_valid]
    
    selected_list_number = (1:length(log_lik_set))[unique_log_lik_order %in% selected_idx]
    
    Q_support_initial_list_list = list()
    
    for (i in 1:length(selected_list_number))
    {
      Q_support_initial_list_list[[i]] = 
        result_per_attempt_whole[[ (selected_list_number[i]) ]]$Q_support_initial_list
    }
    
  return(Q_support_initial_list_list)
}

#Directional Derivative or Gradient function
Directional_mixture = function(phi, x, Q_support_list, gamma_l)
{
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  L = length(Q_support_list)
  
  denominator_density = 0
  for (l_temp in 1:L)
  {
    denominator_density = 
      denominator_density + 
      density_mixture_v2(x = x, 
                         Q_support = Q_support_list[[l_temp]])
  }
  
  pi_fixed = 0
  delta_fixed = 0
  
  for (l_temp in c(1:L)[-gamma_l])
  {
    pi_fixed = pi_fixed + sum(Q_support_list[[l_temp]]$prob)
    delta_fixed = delta_fixed + density_mixture_v2(x = x,
                                                   Q_support = Q_support_list[[l_temp]])
    
  }
  
  if (length(Q_support_list[[gamma_l]]$prob) != 1)
  {
    pi_fixed = pi_fixed + Q_support_list[[gamma_l]]$prob[1]
    delta_fixed = 
      delta_fixed + 
      Q_support_list[[gamma_l]]$prob[1] * 
      dmsn_sahu_v3_biv(x = x,
                       xi = Q_support_list[[gamma_l]]$xi,
                       Sigma = (Q_support_list[[gamma_l]]$alpha[[1]] * 
                                  Q_support_list[[gamma_l]]$Sigma),
                       Lambda = Q_support_list[[gamma_l]]$Lambda)
  }
  
  numerator_density = 
    (1 - pi_fixed)*
    dmsn_sahu_v3_biv(x = x,
                     xi = Q_support_list[[gamma_l]]$xi,
                     Sigma = (phi*
                                Q_support_list[[gamma_l]]$Sigma),
                     Lambda = Q_support_list[[gamma_l]]$Lambda) +
    delta_fixed
  
  result = sum(numerator_density/denominator_density) - n 
  
  return(result)
}

give_interval_with_potential_mixture = function(scale_interval, 
                                                grid_length, 
                                                x, Q_support_list, 
                                                gamma_l, 
                                                maximum_grid_length = 2*grid_length)
{
  grid_for_gradient = seq(from = min(scale_interval), 
                          to = max(scale_interval),
                          length.out = grid_length)
  
  derivatives_of_gradient_function =
    sapply(X = grid_for_gradient, 
           FUN = function(value){
             central_difference(func = Directional_mixture, 
                                value = value, 
                                x = x, gamma_l = gamma_l,
                                Q_support_list = Q_support_list)
           })

  
  derivatives_positive = derivatives_of_gradient_function > 0
  
  grid_length_temp = grid_length
  scale_interval_temp = scale_interval
  
  inc_grid_length_by = ceiling(grid_length/5)
  
  while ( (rev(derivatives_positive)[1] == T) & (length(derivatives_positive) < maximum_grid_length) )
  {
    grid_gap = grid_for_gradient[2] - grid_for_gradient[1]
    
    scale_interval_temp[which.max(scale_interval_temp)] = max(scale_interval_temp) + inc_grid_length_by*grid_gap 
    
    grid_length_old = grid_length_temp
    grid_length_temp = grid_length_temp + inc_grid_length_by
    
    grid_for_gradient = seq(from = min(scale_interval_temp), 
                            to = max(scale_interval_temp),
                            length.out = grid_length_temp) 
    
    additional_grid =
      grid_for_gradient[(grid_length_old+1):(grid_length_old + inc_grid_length_by)]
    
    additional_derivative =
      sapply(X = additional_grid, 
             FUN = function(value){
               central_difference(func = Directional_mixture, 
                                  value = value, 
                                  x = x, gamma_l = gamma_l,
                                  Q_support_list = Q_support_list)
             })
    
    derivatives_positive = c(derivatives_positive, (additional_derivative > 0))
  }
  
  grid_with_potential_mat = where_potential_maximum_grid(logical_vector = derivatives_positive)
  
  result_list = list()
  result_list$grid_for_gradient = grid_for_gradient
  result_list$grid_with_potential_mat = grid_with_potential_mat
  
  return(result_list)
}

CNM_part_mixture = function(x, Q_support_list, gamma_l, 
                            gamma_tol = length(x)/ncol(x) * 1e-14,
                            scale_interval = c(1, 100), 
                            grid_length = 25,
                            tol_Directional_Derivative = 0.01)
{
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  L = length(Q_support_list)
  
  min_scale_interval = min(scale_interval)
  max_scale_interval = max(scale_interval)
  
  max_min_scale_interval = 
    c(max(Q_support_list[[gamma_l]]$alpha), 2*max_scale_interval) %>% min()
  
  max_scale_interval_free = max(c(max_scale_interval, max_min_scale_interval))
  
  scale_interval_free = c(min_scale_interval, max_scale_interval_free)
  
  interval_with_potential =
    give_interval_with_potential_mixture(scale_interval = scale_interval_free,
                                         grid_length = grid_length, 
                                         x = x, Q_support_list = Q_support_list, 
                                         gamma_l = gamma_l)
  
  grid_with_potential_mat = interval_with_potential$grid_with_potential_mat
  grid_for_gradient = interval_with_potential$grid_for_gradient
  
  gradient_results = 
    apply(X = grid_with_potential_mat, MARGIN = 1,
          FUN = function(index){
            optimize(Directional_mixture,
                     interval = c(grid_for_gradient[index[1]], grid_for_gradient[index[2]]),
                     maximum = TRUE, 
                     x = x, Q_support_list = Q_support_list, 
                     gamma_l = gamma_l)
          })
  
  gradient_results_mat = matrix(unlist(gradient_results), ncol = 2, byrow = T)
  colnames(gradient_results_mat) = c("maximum", "objective")
  
  idx_maximizer_among_grids = which.max(gradient_results_mat[,"objective"])
  
  Directional_Derivative_set = gradient_results_mat[, "objective"]
  maximizer_DD_set = gradient_results_mat[, "maximum"]
  
  Directional_Derivative = Directional_Derivative_set[idx_maximizer_among_grids]
  
  if (tol_Directional_Derivative >= Directional_Derivative)
  {
    result_list = list()
    result_list$Q_support_list = Q_support_list
    result_list$Directional_Derivative = Directional_Derivative
    result_list$gamma_l = gamma_l
    return(result_list)
  }
  
  Q_support_list_plus = Q_support_list
  
  prob_Q_support_l = sum(Q_support_list[[gamma_l]]$prob)
  Q_support_list_plus[[gamma_l]]$alpha = c(Q_support_list[[gamma_l]]$alpha, maximizer_DD_set)
  Q_support_list_plus[[gamma_l]]$prob =  c(Q_support_list[[gamma_l]]$prob, rep(0, length(maximizer_DD_set)))
  
  S_l_denominator_vec = 0
  for (l_temp in 1:L)
  {
    S_l_denominator_vec = 
      S_l_denominator_vec + 
      density_mixture_v2(x = x, 
                         Q_support = Q_support_list[[l_temp]])
  }
  
  S_l_denominator = matrix(data = S_l_denominator_vec,
                           nrow = n, ncol = length(Q_support_list_plus[[gamma_l]]$alpha) ) #-1
  
  S_l_numerator = matrix(0, nrow = n, ncol = ncol(S_l_denominator))
  
  for (i in 1:length(Q_support_list_plus[[gamma_l]]$alpha) ) #-1
  {
    S_l_numerator_col_i =
      dmsn_sahu_v3_biv(x = x, 
                       xi = Q_support_list_plus[[gamma_l]]$xi, 
                       Lambda = Q_support_list_plus[[gamma_l]]$Lambda, 
                       Sigma = Q_support_list_plus[[gamma_l]]$alpha[i]*Q_support_list_plus[[gamma_l]]$Sigma)
    
    S_l_numerator[,i] = S_l_numerator_col_i
  }
  
  S_l_mat = S_l_numerator/S_l_denominator
  
  pi_fixed = 0
  for (l_temp in c(1:L)[-gamma_l])
  {
    pi_fixed = pi_fixed + sum(Q_support_list[[l_temp]]$prob)
  }
  
  y_upper = S_l_mat %*% Q_support_list_plus[[gamma_l]]$prob + rep(1, n)
  
  x_4nnls = rbind(S_l_mat*sqrt(gamma_tol), rep(1, ncol(S_l_mat)))
  y_4nnls = c(y_upper*sqrt(gamma_tol), (1-pi_fixed) )
  
  nnls_result = nnls(x_4nnls, y_4nnls)
  
  normalizing_constant = prob_Q_support_l/sum(nnls_result$x)
  Q_support_list_plus[[gamma_l]]$prob = normalizing_constant*nnls_result$x
  
  idx_zero = Q_support_list_plus[[gamma_l]]$prob == 0
  idx_zero[1] = F 
  
  Q_support_list[[gamma_l]]$prob = Q_support_list_plus[[gamma_l]]$prob[!idx_zero]
  Q_support_list[[gamma_l]]$alpha = Q_support_list_plus[[gamma_l]]$alpha[!idx_zero]
  
  result_list = list()
  result_list$Q_support_list = Q_support_list
  result_list$Directional_Derivative = Directional_Derivative
  result_list$gamma_l = gamma_l
  
  return(result_list)
}

EM_part_mixture = function(x, Q_support_list, 
                           tol_EM = 1e-3, max_iter_EM = 100, em_log_lik_check = F)
{
  ##############
  #   MLE(EM)  # Find Fixed Parameters with EM algorithm
  ##############

  L = length(Q_support_list)
  Q_support_list_old = Q_support_list
  Theta_list_old = 
    Theta_list_generator_mixture(Q_support_list = Q_support_list_old)$Theta_list
  
  iter_inner = 0
  log_lik_diff_inner = tol_EM + 1
  log_lik_set = 
    density_mixture_Theta_list(x = x, Theta_list = Theta_list_old) %>% log() %>% sum()
  
  if (em_log_lik_check == T)
  {
    result_prob_em = NULL
    result_xi_em = NULL
    result_Lambda_em = NULL
    result_Sigma_em = NULL
  }
  
  while ( (log_lik_diff_inner > tol_EM) & (iter_inner < max_iter_EM) )
  {
    iter_inner = iter_inner + 1
    
    ##########################
    # MLE(EM) - 1-1. E(Z_ij) #
    ##########################
    ##########################
    # MLE(EM) - 1-2. E(tau) = eta and E(tau %*% t(tau)) = Psi #
    ##########################
    
    Theta_list_result = Theta_list_generator_mixture(Q_support_list = Q_support_list)
    Theta_list = Theta_list_result$Theta_list
    idx_lg_l = Theta_list_result$idx_lg_l
    idx_lg_g = Theta_list_result$idx_lg_g
    number_g = length(idx_lg_g)
    
    #Z
    #tau, eta, Psi 
    Z = Z_generator(x = x, Theta_list = Theta_list)
    mu_tau_TN =  lapply(X = Theta_list, FUN = mu_tau_TN_generator, x = x)
    eta_ij = lapply(X = mu_tau_TN, FUN = eta_generator)
    Psi_ij = lapply(X = mu_tau_TN, FUN = Psi_generator)
    
    ##########################
    # MLE(EM) - 2. M-Step  #
    ##########################
    
    ##########################
    # MLE(EM) - 2-1. Update prob, w  #
    # In addition, put information of Z, eta, Psi
    ##########################
    
    prob_new = apply(Z, 2, FUN = sum)/sum(Z)
    
    Theta_list_new_prob = Theta_list
    for (i in 1:number_g)
    {
      Theta_list_new_prob[[i]]$prob = prob_new[i]
    }
    Q_support_list_new_prob = 
      Q_support_list_generator(Theta_list = Theta_list_new_prob,
                               idx_lg_l = idx_lg_l,
                               idx_lg_g = idx_lg_g)
    Theta_list_new_prob = 
      Theta_list_generator_mixture(Q_support_list = Q_support_list, 
                                   Z = Z, eta_ij = eta_ij, Psi_ij = Psi_ij)$Theta_list
    
    ### Check whether Updating "Prob" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_prob = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_prob) %>% log() %>% sum()
      result_prob_em = c(result_prob_em, Q_result_prob)
    }
    
    ##########################
    # MLE(EM) - 2-2. Update location parameter, xi  #
    ##########################
    
    xi_new = NULL
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      xi_temp =
        xi_generator_v2(x = x,
                        Theta_list = Theta_list_new_prob[idx_lg_only_l])
      
      xi_new = rbind(xi_new, xi_temp)
    }
    
    Theta_list_new_xi = Theta_list_new_prob
    for (i in 1:number_g)
    {
      Theta_list_new_xi[[i]]$xi = xi_new[(idx_lg_l[i]),]
    }
    
    ### Check whether Updating "Xi" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_xi = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_xi) %>% log() %>% sum()
      result_xi_em = c(result_xi_em, Q_result_xi)
    }
    
    ##########################
    # MLE(EM) - 2-3. Update skewness parameter, Lambda#
    ##########################
    Lambda_new = list()
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      Lambda_temp =
        Lambda_generator_v2(x = x, 
                            Theta_list = Theta_list_new_xi[idx_lg_only_l])
      
      Lambda_new[[l_temp]] = Lambda_temp
    }
    
    Theta_list_new_Lambda = Theta_list_new_xi
    for (i in 1:number_g)
    {
      Theta_list_new_Lambda[[i]]$Lambda = Lambda_new[[(idx_lg_l[i])]]
    }
    
    ### Check whether Updating "Lambda" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Lambda = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_Lambda) %>% log() %>% sum()
      result_Lambda_em = c(result_Lambda_em, Q_result_Lambda)
    }
    
    ##########################
    # MLE(EM) - 2-4. Update scale parameter, Sigma#
    ##########################
    
    Sigma_new = list()
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      Sigma_temp = 
        Sigma_generator_v2(x = x,
                           Theta_list = Theta_list_new_Lambda[idx_lg_only_l],
                           Z = Z[, idx_lg_only_l, drop = F])
      
      Sigma_new[[l_temp]] = (Sigma_temp + t(Sigma_temp))/2
    }
    
    Theta_list_new_Sigma = Theta_list_new_Lambda
    for (i in 1:number_g)
    {
      Theta_list_new_Sigma[[i]]$Sigma = 
        Sigma_new[[(idx_lg_l[i])]] * ( Theta_list_new_Sigma[[i]]$alpha_i )
    }
    
    ### Check whether Updating "Sigma" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Sigma = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_Sigma) %>% log() %>% sum()
      result_Sigma_em = c(result_Sigma_em, Q_result_Sigma)
    }
    
    ##########################
    # MLE(EM) - 2-5. Save all updated parameters
    ##########################
    Theta_list = Theta_list_new_Sigma
    Q_support_list = 
      Q_support_list_generator(Theta_list = Theta_list, 
                               idx_lg_l = idx_lg_l, idx_lg_g = idx_lg_g)
    
    ##########################
    # MLE(EM) - 2-6. Calculate to identify wheter stopping rules are satisified or not
    ##########################
    log_lik_new_em = 
      density_mixture_Theta_list(x = x, Theta_list = Theta_list) %>% log() %>% sum()
    
    log_lik_set = c(log_lik_set, log_lik_new_em)
    

    log_lik_diff_inner = log_lik_set[(iter_inner+1)] - log_lik_set[iter_inner]
  }
  
  #Result
  result_list = list()
  
  result_list$Q_support_list = Q_support_list
  
  result_list$iter_inner = iter_inner
  result_list$log_lik_set = log_lik_set
  
  if(em_log_lik_check == T)
  {
    result_list$result_prob_em = result_prob_em
    result_list$result_xi_em = result_xi_em 
    result_list$result_Lambda_em = result_Lambda_em
    result_list$result_Sigma_em = result_Sigma_em
  }
  
  return(result_list)
}

Theta_list_generator_mixture = function(Q_support_list, Z = NULL, eta_ij = NULL, Psi_ij = NULL)
{
  L = length(Q_support_list)
  p = length(Q_support_list[[1]]$xi)
  
  idx_lg_l = NULL
  idx_lg_g = NULL
  alpha_lg = NULL
  
  for (l in 1:L)
  {
    idx_lg_l = c(idx_lg_l, 
                 rep(l, times = length(Q_support_list[[l]]$prob)))
    idx_lg_g = c(idx_lg_g, 1:length(Q_support_list[[l]]$prob))
    
    alpha_lg = c(alpha_lg, Q_support_list[[l]]$alpha)
  }
  
  number_g = length(idx_lg_g)
  
  Theta_list = vector("list", number_g)
  for (i in 1:number_g)
  {
    l = idx_lg_l[i]
    g = idx_lg_g[i]
    
    Theta_list[[i]]$xi = Q_support_list[[l]]$xi
    Theta_list[[i]]$Sigma = Q_support_list[[l]]$Sigma * Q_support_list[[l]]$alpha[g]
    Theta_list[[i]]$Lambda = Q_support_list[[l]]$Lambda
    Theta_list[[i]]$prob = Q_support_list[[l]]$prob[g]
    Theta_list[[i]]$lg = c(l, g)
    
    if (!is.null(Z)) { Theta_list[[i]]$Z_i = Z[, i] }
    if (!is.null(eta_ij)) { Theta_list[[i]]$eta_i = eta_ij[[i]] }
    if (!is.null(Psi_ij)) 
    {
      Theta_list[[i]]$Psi_i =
        Psi_ij[[i]] %>% 
        unlist() %>%
        array(dim = c(p, p, length(Psi_ij[[i]])))
    }
    
    Theta_list[[i]]$alpha_i = Q_support_list[[l]]$alpha[g]
  }
  
  result_list = list()
  result_list$Theta_list = Theta_list
  result_list$idx_lg_l = idx_lg_l
  result_list$idx_lg_g = idx_lg_g
  
  return(result_list)
}

Q_support_list_generator = function(Theta_list, idx_lg_l, idx_lg_g)
{
  Q_support_list = list()
  
  idx_to_be_kept = (1:length(idx_lg_g))[idx_lg_g == 1]
  idx_not_to_be_kept = (1:length(idx_lg_g))[idx_lg_g != 1]
  
  for (l1 in 1:length(idx_to_be_kept))
  {
    Q_support_list[[l1]] = list()
    Q_support_list[[l1]]$alpha = Theta_list[[(idx_to_be_kept[l1])]]$alpha_i
    Q_support_list[[l1]]$xi = Theta_list[[(idx_to_be_kept[l1])]]$xi
    Q_support_list[[l1]]$Sigma = 
      Theta_list[[(idx_to_be_kept[l1])]]$Sigma / Theta_list[[(idx_to_be_kept[l1])]]$alpha_i
    Q_support_list[[l1]]$Lambda = Theta_list[[(idx_to_be_kept[l1])]]$Lambda
    Q_support_list[[l1]]$prob = Theta_list[[(idx_to_be_kept[l1])]]$prob
  }
  
  for (l_not1 in idx_not_to_be_kept)
  {
    l = idx_lg_l[ l_not1 ]
    g = idx_lg_g[ l_not1 ]
    
    Q_support_list[[l]]$alpha = 
      c(Q_support_list[[l]]$alpha, Theta_list[[l_not1]]$alpha_i)
    Q_support_list[[l]]$prob = 
      c(Q_support_list[[l]]$prob, Theta_list[[l_not1]]$prob)
  }
  
  return(Q_support_list)
}

density_mixture_Theta_list = function(x, Theta_list)
{
  number_g = length(Theta_list)
  n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  
  density_list = 
    lapply(X = Theta_list,
           FUN = function(Theta_list_i){ # correct
             i_density = 
               (Theta_list_i$prob*
                  dmsn_sahu_v3_biv(x = x, 
                                   xi = Theta_list_i$xi,
                                   Sigma = Theta_list_i$Sigma, Lambda = Theta_list_i$Lambda))
           })
  
  density_vector = 
    matrix(as.numeric(unlist(density_list)), nrow = n, ncol = number_g) %>% 
    apply(MARGIN = 1, FUN = sum)
  
  return(density_vector)
}

density_mixture_Q_support_list = function(x, Q_support_list)
{
  lik_list = 
    lapply(X = Q_support_list,
           FUN = function(Q_support_i){
             density_mixture_v2(x = x, Q_support = Q_support_i)
           })
  
  lik_mat = lik_list %>% unlist %>% matrix(ncol = length(Q_support_list))
  lik_vec = lik_mat %>% apply(MARGIN = 1, FUN = sum)
  
  return(lik_vec)
}

Z_giver = function(x, Q_support_list)
{
  Theta_list_result = Theta_list_generator_mixture(Q_support_list = Q_support_list)
  
  Theta_list = Theta_list_result$Theta_list
  idx_lg_l = Theta_list_result$idx_lg_l
  
  L = length(Q_support_list)
  number_g = length(Theta_list)
    
  Z = Z_generator(x = x, Theta_list = Theta_list)
  
  Z_outer = matrix(0, nrow = nrow(Z), ncol = L)
  
  for (l_temp in 1:L)
  {
    idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
    
    Z_l = Z[,idx_lg_only_l, drop = F]
    
    Z_outer[,l_temp] = Z_l %>% apply(1, sum)
  }
  
  return(Z_outer)
}


CNM_SSMMSN_mixture = function(X, Q_support_initial_list,
                              tol_outer, tol_Directional_Derivative = 0.1,
                              tol_EM = tol_outer,
                              max_iter_outer = 10, max_iter_EM = 1, max_iter_CNM = 1,
                              scale_interval = c(1, 100), grid_length = 100,
                              print_each_iteration = T)
{
  Q_support_list = Q_support_initial_list
  
  log_lik_diff_outer = tol_outer + 1
  Directional_Derivative = tol_Directional_Derivative + 1
  iter = 0
  
  LL_set = density_mixture_Q_support_list(x = X, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
  DD_mat = NULL
  
  while((iter < max_iter_outer) & (log_lik_diff_outer > tol_outer)) # 
  {
    start_time = Sys.time()
    iter = iter + 1
    
    Log_Likelihood_old = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    # EM
    EM_result = 
      EM_part_mixture_decme(x = X, Q_support_list = Q_support_list, 
                      max_iter_EM = max_iter_EM, tol_EM = tol_EM, em_log_lik_check = F)
    
    Q_support_list = EM_result$Q_support_list
    
    L_inter = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    Directional_set = NULL
    Support_set = NULL
    iter_CNM_set = NULL
    # CNM
    for (l in 1:length(Q_support_list))
    {
      iter_CNM = 0
      iter_switch = T
      while ((iter_CNM < max_iter_CNM) & (iter_switch))
      {
        iter_CNM = iter_CNM + 1
        CNM_result = 
          CNM_part_mixture(x = X, Q_support_list = Q_support_list, 
                           gamma_l = l,
                           scale_interval = scale_interval, grid_length = grid_length, 
                           tol_Directional_Derivative = tol_Directional_Derivative)
        
        if (CNM_result$Directional_Derivative >= tol_Directional_Derivative)
        {
          Q_support_list = CNM_result$Q_support_list
        } else {
          iter_switch = F
          next 
        }
      }
      
      iter_CNM_set = c(iter_CNM_set, iter_CNM)
      Directional_set = c(Directional_set, CNM_result$Directional_Derivative)
      Support_set = c(Support_set, length(Q_support_list[[l]]$alpha) )
    }
    
    Log_Likelihood = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    spent_time = Sys.time() - start_time
    
    LL_set = c(LL_set, Log_Likelihood)
    DD_mat = rbind(DD_mat, Directional_set)
    
    log_lik_diff_outer = LL_set[iter + 1] - LL_set[iter]
    
    if (print_each_iteration == T)
    {
      cat("########################################################", "\n")
      cat("iteration spent for EM = ", EM_result$iter_inner, "\n")
      cat("iteration spent for CNM = ", paste(iter_CNM_set, collapse = ", "), "\n")
      cat("iteration = " , iter , 
          "D_set =", paste(round(Directional_set, 2), collapse = ", "),  
          "L = ", Log_Likelihood, 
          "L_inter = ", L_inter,
          "L_old = ", Log_Likelihood_old, 
          "Support =", paste(Support_set, collapse = ", "), "\n")
      cat("Spent time: ") 
      print(spent_time) 
    }
  }
  
  result_list = list()
  result_list$Q_support_list = Q_support_list
  result_list$LL_set = LL_set
  result_list$DD_mat = DD_mat
  
  return(result_list)
}

#CNM_SMSN
CNM_SMSN_mixture = function(X, Q_support_initial_list,
                            tol_outer = 1e-4, tol_Directional_Derivative = 1e-4,
                            tol_EM = tol_outer,
                            max_iter_outer = 100, max_iter_EM = 5, max_iter_CNM = 5, 
                            scale_interval = c(1, 100), grid_length = 25,
                            print_each_iteration = F)
{
  Q_support_list = Q_support_initial_list
  
  log_lik_diff_outer = tol_outer + 1
  Directional_Derivative = tol_Directional_Derivative + 1
  iter = 0
  
  if (print_each_iteration == T) 
  {
    print(Sys.time())
    print("Biscuit firing on-going")
  }
  
  bisque_firing_result =
    EM_part_mixture_decme(x = X, Q_support_list = Q_support_list, 
                          max_iter_EM = max_iter_outer*0.8, 
                          tol_EM = (tol_EM*5), em_log_lik_check = F)
  
  if (print_each_iteration == T) 
  {
    print(Sys.time())
    cat("Biscuit firing is done with", bisque_firing_result$iter_inner, "iterations \n")
  }
  
  LL_set = bisque_firing_result$log_lik_set
  DD_mat = NULL
  iter = bisque_firing_result$iter_inner
  Q_support_list = bisque_firing_result$Q_support_list
  
  while((iter < max_iter_outer) & (log_lik_diff_outer > tol_outer)) # 
  {
    start_time = Sys.time()
    iter = iter + 1
    
    Log_Likelihood_old = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    # CNM
    Directional_set = NULL
    Support_set = NULL
    iter_CNM_set = NULL
    for (l in 1:length(Q_support_list))
    {
      iter_CNM = 0
      iter_switch = T
      while ((iter_CNM < max_iter_CNM) & (iter_switch))
      {
        iter_CNM = iter_CNM + 1
        CNM_result = 
          CNM_part_mixture(x = X, Q_support_list = Q_support_list, 
                           gamma_l = l,
                           scale_interval = scale_interval, grid_length = grid_length, 
                           tol_Directional_Derivative = tol_Directional_Derivative)
        
        if (CNM_result$Directional_Derivative >= tol_Directional_Derivative)
        {
          Q_support_list = CNM_result$Q_support_list
        } else {
          iter_switch = F
          next 
        }
      }
      
      iter_CNM_set = c(iter_CNM_set, iter_CNM)
      Directional_set = c(Directional_set, CNM_result$Directional_Derivative)
      Support_set = c(Support_set, length(Q_support_list[[l]]$alpha) )
    }
    
    L_inter = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    # EM
    EM_result = 
      EM_part_mixture_decme(x = X, Q_support_list = Q_support_list, 
                            max_iter_EM = max_iter_EM, tol_EM = tol_EM, em_log_lik_check = F)
    
    Q_support_list = EM_result$Q_support_list
    
    Log_Likelihood = 
      density_mixture_Q_support_list(x = X, Q_support_list = Q_support_list) %>% log() %>% sum()
    
    spent_time = Sys.time() - start_time
    
    LL_set = c(LL_set, Log_Likelihood)
    DD_mat = rbind(DD_mat, Directional_set)

    log_lik_diff_outer = LL_set[iter + 1] - LL_set[iter]
    
    if (print_each_iteration == T)
    {
      cat("########################################################", "\n")
      cat("iteration spent for EM = ", EM_result$iter_inner, "\n")
      cat("iteration spent for CNM = ", paste(iter_CNM_set, collapse = ", "), "\n")
      cat("iteration = " , iter , 
          "D_set =", paste(round(Directional_set, 2), collapse = ", "),  
          "L = ", Log_Likelihood, 
          "L_inter = ", L_inter,
          "L_old = ", Log_Likelihood_old, 
          "Support =", paste(Support_set, collapse = ", "), "\n")
      cat("Spent time: ") 
      print(spent_time) 
    }
  }
  
  result_list = list()
  result_list$Q_support_list = Q_support_list
  result_list$LL_set = LL_set
  result_list$DD_mat = DD_mat
  
  return(result_list)
}

### Simulation Code

simulation_function_several_initial = function(X, L, log_lik_tol, max_iter, label_true,
                                               initial_label_count_minimum = 5,
                                               how_much_initial = 5,
                                               zero_Lambda = T)
{
  tryCatch(
    expr =
      {
        if (zero_Lambda == T)
        {
          how_much_initial = how_much_initial*2
          Q_support_initial_list_list =
            FMSN_initial_with_zero_Lambda_5clustering(x = X, L = L, 
                                                      how_much_initial = how_much_initial, how_much_per_clustering = 1, 
                                                      label_count_minimum = initial_label_count_minimum)
        } else if (zero_Lambda == F){
          Q_support_initial_list_list =
            FMSN_initial_5clustering(x = X, L = L, 
                                     how_much_initial = how_much_initial, how_much_per_clustering = 1, 
                                     label_count_minimum = initial_label_count_minimum)
        } else {
          print("Wrong zero lambda option")
          return("wrong")
        }
        
        
        ## 1. F - SMSN
        
        result_FSMSN_list =
          lapply(X = Q_support_initial_list_list,
                 FUN = function(Q_support_initial_list){
                   result_FSMSN = list()
                   result_FSMSN$LL_set = -Inf
                   start_time_FSMSN = Sys.time()
                   
                   tryCatch(expr = {
                     result_FSMSN =
                       CNM_SSMMSN_mixture_V2(X = X, Q_support_initial_list = Q_support_initial_list,
                                             tol_outer = log_lik_tol, 
                                             tol_Directional_Derivative = 1e-2, # 1e-2
                                             max_iter_outer = max_iter, max_iter_EM = 5, max_iter_CNM = 5,
                                             scale_interval = c(1, 100), grid_length = 25,
                                             print_each_iteration = F)
                   }, error = function(e){print(e)},
                   warning = function(w) {print(w)})
                   
                   spent_time_FSMSN = Sys.time() - start_time_FSMSN
                   result_FSMSN$spent_time = spent_time_FSMSN
                   
                   return(result_FSMSN)
                 }
          )
        
        
        log_lik_FSMSN = 
          lapply(X = result_FSMSN_list,
                 FUN = function(result_FSMSN){
                   return(rev(result_FSMSN$LL_set)[1])
                 }) %>% unlist()
        
        result_FSMSN = result_FSMSN_list[[ (log_lik_FSMSN %>% which.max()) ]]
        spent_time_FSMSN = result_FSMSN$spent_time
        
        X_Z_FSMSN = Z_giver(x = X, Q_support_list = result_FSMSN$Q_support_list)
        label_est_FSMSN = apply(X_Z_FSMSN, 1, which.max) %>% as.character()
        
        log_lik = rev(result_FSMSN$LL_set)[1]
        p = result_FSMSN$Q_support_list[[1]]$xi %>% length()
        n = nrow(X)
        num_para = 
          (unlist(result_FSMSN$Q_support_list) %>% length()) - L - (L*p*(p-1)) 
        # scale parameter들 하나씩은 fixed니까 L개 만큼 빼준다
        # 그리고 Lambda에 대각행렬 제약조건 있으니 반영해서 빼준다
        # proportion도 parameter 니까 넣어주려고 했는데 list에 이미 들어간 정보구나 
        
        BIC_FSMSN = -2*log_lik + (num_para*log(n))
        
        # 2. F - MSN
        result_FMSN_list =
          lapply(X = Q_support_initial_list_list,
                 FUN = function(Q_support_initial_list){
                   result_FMSN = list()
                   result_FMSN$log_lik_set = -Inf
                   start_time_FMSN = Sys.time()
                   
                   tryCatch(expr = {
                   result_FMSN = 
                     EM_part_mixture_decme(x = X, Q_support_list = Q_support_initial_list,
                                     tol_EM = log_lik_tol, max_iter_EM = max_iter, 
                                     em_log_lik_check = F)
                   }, error = function(e){print(e)},
                   warning = function(w) {print(w)})
                   
                   spent_time_FMSN = Sys.time() - start_time_FMSN
                   result_FMSN$spent_time = spent_time_FMSN
                   
                   return(result_FMSN)
                 })
        
        log_lik_FMSN = 
          lapply(X = result_FMSN_list,
                 FUN = function(result_FMSN){
                   return(result_FMSN$log_lik_set %>% max())
                 }) %>% unlist()
        
        result_FMSN = result_FMSN_list[[ (log_lik_FMSN %>% which.max()) ]]
        spent_time_FMSN = result_FMSN$spent_time
        
        X_Z_FMSN = Z_giver(x = X, Q_support_list = result_FMSN$Q_support_list)
        label_est_FMSN = apply(X_Z_FMSN, 1, which.max) %>% as.character()
        
        log_lik = rev(result_FMSN$log_lik_set)[1]
        num_para = (p + p^2 + p + 1)*L #proportion 도 파라미터 이기 때문이다.
        
        BIC_FMSN = -2*log_lik + (num_para*log(n))
        
      
        ## 3. F - MVN
        start_time_FMVN = Sys.time()
        suppressPackageStartupMessages(expr = library(mclust))
        result_FMVN = mclust::Mclust(data = X, G = L, control = emControl(tol = log_lik_tol), verbose = F)
        detach("package:mclust", unload=TRUE)
        spent_time_FMVN = Sys.time() - start_time_FMVN
        
        label_est_FMVN = result_FMVN$classification
        
        # -2*result_FMVN$loglik + (6*L - 6)*log(n) #6을 왜 빼지
        
        BIC_FMVN = -result_FMVN$bic
        
        ## 4. kmeans
        start_time_kmeans = Sys.time()
        result_kmeans = kmeans(x = X, centers = L, nstart = 40)
        spent_time_kmeans =  Sys.time() - start_time_kmeans
        
        label_est_kmeans = result_kmeans$cluster 
        
        BIC_kmeans = NA
        
        result_list = list()
        
        # result_list
        result_list$result_list = list()
        result_list$result_list$result_FSMSN = result_FSMSN
        result_list$result_list$result_FMSN = result_FMSN
        result_list$result_list$result_FMVN = result_FMVN
        result_list$result_list$result_kmeans = result_kmeans
        
        # ARI, AMI
        ARI_FSMSN = ARI(c1 = label_est_FSMSN, c2 = label_true)
        ARI_FMSN = ARI(c1 = label_est_FMSN, c2 = label_true)
        ARI_FMVN = ARI(c1 = label_est_FMVN, c2 = label_true)
        ARI_kmeans = ARI(c1 = label_est_kmeans, c2 = label_true)
        
        AMI_FSMSN = AMI(c1 = label_est_FSMSN, c2 = label_true)
        AMI_FMSN = AMI(c1 = label_est_FMSN, c2 = label_true)
        AMI_FMVN = AMI(c1 = label_est_FMVN, c2 = label_true)
        AMI_kmeans = AMI(c1 = label_est_kmeans, c2 = label_true)
        
        ARI_AMI_BIC_time_mat =
          matrix(data =
                   c(ARI_FSMSN, ARI_FMSN, ARI_FMVN, ARI_kmeans,
                     AMI_FSMSN, AMI_FMSN, AMI_FMVN, AMI_kmeans,
                     BIC_FSMSN, BIC_FMSN, BIC_FMVN, BIC_kmeans,
                     spent_time_FSMSN, spent_time_FMSN, spent_time_FMVN, spent_time_kmeans),
                 byrow = T, nrow = 4)
        
        colnames(ARI_AMI_BIC_time_mat) = c("FSMSN", "FMSN", "FMVN", "kmeans")
        rownames(ARI_AMI_BIC_time_mat) = c("ARI", "AMI", "BIC", "Time")
        
        result_list$ARI_AMI_BIC_time_mat = ARI_AMI_BIC_time_mat
        
        result_list$X = X
        
        return(result_list)
      }, 
    error = function(e){print(e)}
  )
}

# 
EM_part_mixture_decme = function(x, Q_support_list, 
                                 tol_EM = 1e-3, 
                                 max_iter_EM = 100, em_log_lik_check = F,
                                 interval_line_search = c(-10, 10))
{
  ##############
  #   MLE(EM)  # Find Fixed Parameters with EM algorithm
  ##############
  
  #n = ifelse(is.null(dim(x)), yes = length(x), no = nrow(x))
  L = length(Q_support_list)
  Q_support_list_old = Q_support_list
  Theta_list_old = 
    Theta_list_generator_mixture(Q_support_list = Q_support_list_old)$Theta_list
  
  iter_inner = 0
  log_lik_diff_inner = tol_EM + 1
  log_lik_set = 
    density_mixture_Theta_list(x = x, Theta_list = Theta_list_old) %>% log() %>% sum()
  
  if (em_log_lik_check == T)
  {
    result_prob_em = NULL
    result_xi_em = NULL
    result_Lambda_em = NULL
    result_Sigma_em = NULL
  }
  #DECME
  Lambda_history = list()
  
  for (l in 1:L)
  {
    Lambda_zero = Q_support_list_old[[l]]$Lambda %>% diag()
    Lambda_history[[l]] = matrix(Lambda_zero, nrow = 1)
  }
  
  
  while ( (log_lik_diff_inner > tol_EM) & (iter_inner < max_iter_EM) )
  {
    iter_inner = iter_inner + 1
    
    ##########################
    # MLE(EM) - 1-1. E(Z_ij) #
    ##########################
    ### Z_ij, Remember that index j is indicating j_th observation.
    ##########################
    # MLE(EM) - 1-2. E(tau) = eta and E(tau %*% t(tau)) = Psi #
    ##########################
    
    Theta_list_result = Theta_list_generator_mixture(Q_support_list = Q_support_list)
    Theta_list = Theta_list_result$Theta_list
    idx_lg_l = Theta_list_result$idx_lg_l
    idx_lg_g = Theta_list_result$idx_lg_g
    number_g = length(idx_lg_g)
    
    #Z
    #tau, eta, Psi 
    Z = Z_generator(x = x, Theta_list = Theta_list)
    mu_tau_TN =  lapply(X = Theta_list, FUN = mu_tau_TN_generator, x = x)
    eta_ij = lapply(X = mu_tau_TN, FUN = eta_generator)
    Psi_ij = lapply(X = mu_tau_TN, FUN = Psi_generator)
    
    ##########################
    # MLE(EM) - 2. M-Step  #
    ##########################
    
    ##########################
    # MLE(EM) - 2-1. Update prob, w  #
    # In addition, put information of Z, eta, Psi
    ##########################
    
    prob_new = apply(Z, 2, FUN = sum)/sum(Z)
    
    Theta_list_new_prob = Theta_list
    for (i in 1:number_g)
    {
      Theta_list_new_prob[[i]]$prob = prob_new[i]
    }
    Q_support_list_new_prob = 
      Q_support_list_generator(Theta_list = Theta_list_new_prob,
                               idx_lg_l = idx_lg_l,
                               idx_lg_g = idx_lg_g)
    Theta_list_new_prob = 
      Theta_list_generator_mixture(Q_support_list = Q_support_list, 
                                   Z = Z, eta_ij = eta_ij, Psi_ij = Psi_ij)$Theta_list
    
    ### Check whether Updating "Prob" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_prob = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_prob) %>% log() %>% sum()
      result_prob_em = c(result_prob_em, Q_result_prob)
    }
    
    ##########################
    # MLE(EM) - 2-2. Update location parameter, xi  #
    ##########################
    
    xi_new = NULL
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      xi_temp =
        xi_generator_v2(x = x,
                        Theta_list = Theta_list_new_prob[idx_lg_only_l])
      
      xi_new = rbind(xi_new, xi_temp)
    }
    
    Theta_list_new_xi = Theta_list_new_prob
    for (i in 1:number_g)
    {
      Theta_list_new_xi[[i]]$xi = xi_new[(idx_lg_l[i]),]
    }
    
    ### Check whether Updating "Xi" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_xi = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_xi) %>% log() %>% sum()
      result_xi_em = c(result_xi_em, Q_result_xi)
    }
    
    ##########################
    # MLE(EM) - 2-3. Update scale parameter, Sigma# This is used to be 2-4
    ##########################
    
    Sigma_new = list()
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      Sigma_temp = 
        Sigma_generator_v2(x = x,
                           Theta_list = Theta_list_new_xi[idx_lg_only_l],
                           Z = Z[, idx_lg_only_l, drop = F])
      
      Sigma_new[[l_temp]] = (Sigma_temp + t(Sigma_temp))/2
    }
    
    Theta_list_new_Sigma = Theta_list_new_xi
    for (i in 1:number_g)
    {
      Theta_list_new_Sigma[[i]]$Sigma = 
        Sigma_new[[(idx_lg_l[i])]] * ( Theta_list_new_Sigma[[i]]$alpha_i )
    }
    
    ### Check whether Updating "Sigma" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Sigma = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_Sigma) %>% log() %>% sum()
      result_Sigma_em = c(result_Sigma_em, Q_result_Sigma)
    }
    
    
    ########################## DECME
    # MLE(EM) - 2-4. Update skewness parameter, Lambda# This is used to be 2-3
    ##########################
    Q_support_list_new_Sigma = Q_support_list_generator(Theta_list = Theta_list_new_Sigma,
                                                        idx_lg_l = idx_lg_l, idx_lg_g = idx_lg_g)
    
    Q_support_list_SOR = Q_support_list_new_Sigma
    
    Lambda_new = list()
    for (l_temp in 1:L)
    {
      idx_lg_only_l = (1:number_g)[idx_lg_l == l_temp]
      
      Lambda_temp =
        Lambda_generator_v2(x = x, 
                            Theta_list = Theta_list_new_Sigma[idx_lg_only_l])
      
      Lambda_history_mat = Lambda_history[[l_temp]]
      
      Lambda_t_minus = Lambda_history_mat[nrow(Lambda_history_mat), ] %>% diag()
      if (nrow(Lambda_history_mat) >= 2) {
        Lambda_t_minus2 = Lambda_history_mat[(nrow(Lambda_history_mat) - 1), ] %>% diag()
      } else {
        Lambda_t_minus2 = Lambda_history_mat[nrow(Lambda_history_mat), ] %>% diag()
      }
      
      d_t_1 = Lambda_temp - Lambda_t_minus
      
      D_CM_step_1 = 
        optimize(f = function(a)
        {
          Lambda_SOR_temp = Lambda_temp + a*d_t_1
          
          Q_support_list_SOR[[l_temp]]$Lambda = Lambda_SOR_temp
          
          log_lik_SOR = 
            density_mixture_Q_support_list(x = x, Q_support_list = Q_support_list_SOR) %>% 
            log() %>% sum()
          
          return(log_lik_SOR)
        }, maximum = T, interval = interval_line_search)
      
      a_t_1 = D_CM_step_1$maximum
      
      Lambda_SOR = Lambda_temp + a_t_1*d_t_1
      
      d_t_2 = Lambda_SOR - Lambda_t_minus2
      
      D_CM_step_2 =
        optimize(f = function(a)
        {
          Lambda_tilde_temp = Lambda_SOR + a*d_t_2 
          
          Q_support_list_SOR[[l_temp]]$Lambda = Lambda_tilde_temp
          
          log_lik_tilde = 
            density_mixture_Q_support_list(x = x, Q_support_list = Q_support_list_SOR) %>% 
            log() %>% sum()
          
          return(log_lik_tilde)
        }, maximum = T, interval = interval_line_search)
      
      a_t_2 = D_CM_step_2$maximum
      
      Lambda_tilde = Lambda_SOR + a_t_2*d_t_2 
      
      Q_support_list_SOR[[l_temp]]$Lambda = Lambda_tilde
      
      Lambda_history[[l_temp]] = rbind(Lambda_history[[l_temp]], 
                                       diag(Lambda_tilde))
      Lambda_new[[l_temp]] = Lambda_tilde
    }
    
    Theta_list_new_Lambda = Theta_list_new_Sigma
    for (i in 1:number_g)
    {
      Theta_list_new_Lambda[[i]]$Lambda = Lambda_new[[(idx_lg_l[i])]]
    }
    
    ### Check whether Updating "Lambda" increases Log-Lik
    if (em_log_lik_check == T)
    {
      Q_result_Lambda = 
        density_mixture_Theta_list(x = x, Theta_list = Theta_list_new_Lambda) %>% log() %>% sum()
      result_Lambda_em = c(result_Lambda_em, Q_result_Lambda)
    }
    
    
    ##########################
    # MLE(EM) - 2-5. Save all updated parameters
    ##########################
    Theta_list = Theta_list_new_Lambda
    Q_support_list = 
      Q_support_list_generator(Theta_list = Theta_list, 
                               idx_lg_l = idx_lg_l, idx_lg_g = idx_lg_g)
    
    ##########################
    # MLE(EM) - 2-6. Calculate to identify wheter stopping rules are satisified or not
    ##########################
    log_lik_new_em = 
      density_mixture_Theta_list(x = x, Theta_list = Theta_list) %>% log() %>% sum()
    
    log_lik_set = c(log_lik_set, log_lik_new_em)
    
    log_lik_diff_inner = log_lik_set[(iter_inner+1)] - log_lik_set[iter_inner] 
  }
  
  #Result
  result_list = list()
  
  result_list$Q_support_list = Q_support_list
  
  result_list$iter_inner = iter_inner
  result_list$log_lik_set = log_lik_set
  
  if(em_log_lik_check == T)
  {
    result_list$result_prob_em = result_prob_em
    result_list$result_xi_em = result_xi_em 
    result_list$result_Sigma_em = result_Sigma_em
    result_list$result_Lambda_em = result_Lambda_em
  }
  
  return(result_list)
}


fmmst_to_Q_support_initial_list_list = function(fmmst_initial_list_list)
{
  Q_support_initial_list_list = 
    lapply(X = fmmst_initial_list_list,
           FUN = function(fmmst_initial_list){
             Q_support_list = vector("list", length(fmmst_initial_list$mu))
             
             for (k in 1:length(fmmst_initial_list$mu))
             {
               Q_support_list[[k]]$alpha = 1
               Q_support_list[[k]]$xi = fmmst_initial_list$mu[[k]] %>% as.vector()
               Q_support_list[[k]]$Sigma = fmmst_initial_list$sigma[[k]]
               Q_support_list[[k]]$Lambda = fmmst_initial_list$delta[[k]] %>% as.vector() %>% diag()
               Q_support_list[[k]]$prob = fmmst_initial_list$pro[k]
             }
             
             return(Q_support_list)
           })
  
  return(Q_support_initial_list_list)
}

#FMSN_initial_generator_true
Q_support_initial_generator = function(x, label_est)
{
  Q_support_initial_list = list()
  
  unique_label = unique(label_est)
  
  for(l in 1:length(unique_label))
  {
    x_l_temp = x[label_est == unique_label[l], , drop = F]
    
    x_l_initial = MSN_initial_generator(x = x_l_temp)
    
    Q_support_l_initial = list()
    Q_support_l_initial$alpha = c(1)
    Q_support_l_initial$xi = x_l_initial$xi_init
    Q_support_l_initial$Sigma = x_l_initial$Sigma_init
    Q_support_l_initial$Lambda = x_l_initial$Lambda_init
    Q_support_l_initial$prob = (label_est == unique_label[l]) %>% mean()
    
    Q_support_initial_list[[l]] = Q_support_l_initial
  }
  
  log_lik_temp = density_mixture_Q_support_list(x = x, Q_support_list = Q_support_initial_list) %>% log() %>% sum()
  
  return(Q_support_initial_list)
}

