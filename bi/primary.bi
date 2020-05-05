model primary {
  const N_1 = 1045
  const N_2 = 2666
  const c_11 = 1
  const X = 0.1 
  const gamma_a = 1 / 5
  const gamma_p = 1 / 2.1
  const gamma_s = 1 / 2.9
  const mu = 1
  const eta = 1 / 7
  const t_mu = 16
  const phi = 199 / 314
  const alpha = 1
  const nu = 2 * (1 / 4.3)
  const b_2 = 10

  param beta_bar
  param c_22
  param b_1
  param tau_11
  param tau_22
  param chi
  param theta_a
  param theta_p

  input dN_tests

  state S_1     (has_output = 0) 
  state S_2     (has_output = 0) 
  state E_1     (has_output = 0)   
  state E_2     (has_output = 0) 
  state E_1e    (has_output = 0) 
  state E_2e    (has_output = 0) 
  state I_a1    (has_output = 0) 
  state I_a2    (has_output = 0)
  state I_p1    (has_output = 0) 
  state I_p2    (has_output = 0) 
  state I_su1   (has_output = 0) 
  state I_su2   (has_output = 0) 
  state I_sk1   (has_output = 0) 
  state I_sk2   (has_output = 0) 
  state T_1     (has_output = 0) 
  state T_2     (has_output = 0) 
  state C_1     (has_output = 0) 
  state C_2     (has_output = 0) 
  
  state Z_sk1   (has_output = 0) 
  state Z_sk2   (has_output = 0) 
  state Z_prev (has_output = 0) 
  
  state f           (has_output = 0)
  state delta_I_a1  (has_output = 0) 
  state delta_I_a2  (has_output = 0) 
  state delta_I_p1  (has_output = 0) 
  state delta_I_p2  (has_output = 0) 
  state delta_T_1   (has_output = 0 
  state delta_T_2   (has_output = 0) 
  
  obs symp_1
  obs symp_2
  obs non_symp
  
  sub initial {
    S_1 <- N_1
    S_2 <- N_2-1
    E_1 <- 0
    E_2 <- 0
    E_1e <- 0
    E_2e <- 0
    I_a1 <- 0
    I_a2 <- 0
    I_p1 <- 0
    I_p2 <- 0
    I_su1 <- 0
    I_su2 <- 0
    I_sk1 <- 0
    I_sk2 <- 1
    T_1 <- 0
    T_2 <- 0
    C_1 <- 0
    C_2 <- 0
    
    Z_sk1 <- 0
    Z_sk2 <- 0
    Z_prev <- 0
  }

  sub parameter {
    beta_bar ~ uniform(0,100)
    c_22 ~ uniform(0,100)
    b_1 ~ uniform(0, 1)
    tau_11 ~ uniform(0, 32)
    tau_22 ~ uniform(0, 32)
    chi ~ uniform(0, 1)
    theta_a ~ uniform(0, 1)
    theta_p ~ uniform(0, 1)
  }

  sub transition (delta = 1.0) {

    inline beta_t11 = beta_bar*(1-b_1/(1+exp(-b_2*(t_now-tau_11))))
    inline beta_t12 = beta_bar*(1-b_1/(1+exp(-b_2*(t_now-tau_22))))
    inline beta_t22 = beta_bar*(1-b_1/(1+exp(-b_2*(t_now-tau_22))))
    
    inline mu_t = (t_now < t_mu) ? 0 : mu
    inline gamma_s_t = (t_now < t_mu) ? gamma_s : 0
    
    inline c_12 = X * c_22

    Z_sk1 <- (t_now % 1 == 0 ? 0 : Z_sk1)
    Z_sk2 <- (t_now % 1 == 0 ? 0 : Z_sk2)
    
    ode {
      dS_1/dt = -beta_t11 * c_11 * ((theta_a * I_a1 + theta_p * I_p1 + I_sk1 + I_su1) / N_1) * S_1 - beta_t12 * c_12 * ((theta_a * I_a2 + theta_p * I_p2 + I_sk2 + I_su2) / N_2) * S_1
      
      dS_2/dt = -beta_t12 * c_12 * ((theta_a * I_a1 + theta_p * I_p1 + I_sk1 + I_su1) / N_1) * S_2 - beta_t22 * c_22 * ((theta_a * I_a2 + theta_p * I_p2 + I_sk2 + I_su2) / N_2) * S_2  
      
      dE_1/dt = beta_t11 * c_11 * ((theta_a * I_a1 + theta_p * I_p1 + I_sk1 + I_su1) / N_1) * S_1 + beta_t12 * c_12 * ((theta_a * I_a2 + theta_p * I_p2 + I_sk2 + I_su2) / N_2) * S_1 - nu * E_1
      dE_2/dt = beta_t12 * c_12 * ((theta_a * I_a1 + theta_p * I_p1 + I_sk1 + I_su1) / N_1) * S_2 + beta_t22 * c_22 * ((theta_a * I_a2 + theta_p * I_p2 + I_sk2 + I_su2) / N_2) * S_2 - nu * E_2
      
      dE_1e/dt = nu * E_1 - nu * E_1e
      dE_2e/dt = nu * E_2 - nu * E_2e
      
      dI_a1/dt = chi * nu * E_1e - gamma_a * I_a1
      dI_a2/dt = chi * nu * E_2e - gamma_a * I_a2
      
      dI_p1/dt = (1 - chi) * nu * E_1e - gamma_p * I_p1
      dI_p2/dt = (1 - chi) * nu * E_2e - gamma_p * I_p2
      
      dI_su1/dt = (1 - phi) * gamma_p * I_p1 - gamma_s_t * I_su1 - mu_t * I_su1
      dI_su2/dt = (1 - phi) * gamma_p * I_p2 - gamma_s_t * I_su2 - mu_t * I_su2
      
      dI_sk1/dt = phi * gamma_p * I_p1 - gamma_s_t * I_sk1 - mu_t * I_sk1
      dI_sk2/dt = phi * gamma_p * I_p2 - gamma_s_t * I_sk2 - mu_t * I_sk2
    
      dT_1/dt = alpha * gamma_a * I_a1 + alpha * gamma_s_t * I_sk1 + alpha * gamma_s_t * I_su1 - eta * T_1
      dT_2/dt = alpha * gamma_a * I_a2 + alpha * gamma_s_t * I_sk2 + alpha * gamma_s_t * I_su2 - eta * T_2
      
      dC_1/dt = (1 - alpha) * gamma_a * I_a1 + eta * T_1 + (1 - alpha) * gamma_s_t * I_sk1 + (1 - alpha) * gamma_s_t * I_su1
      dC_2/dt = (1 - alpha) * gamma_a * I_a2 + eta * T_2 + (1 - alpha) * gamma_s_t * I_sk2 + (1 - alpha) * gamma_s_t * I_su2
      
      dZ_sk1/dt = phi * gamma_p * I_p1
      dZ_sk2/dt = phi * gamma_p * I_p2
    }
    
    Z_prev <- (I_a1 + I_a2 + I_p1 + I_p2 + T_1 + T_2) / (S_1 + S_2 + E_1 + E_2 + E_1e + E_2e + I_a1 + I_a2 + I_p1 + I_p2 + T_1 + T_2 + C_1 + C_2)
    
    f <- dN_tests / (S_1 + S_2 + E_1 + E_2 + E_1e + E_2e + I_a1 + I_a2 + I_p1 + I_p2 + T_1 + T_2 + C_1 + C_2)
  
    delta_I_a1 <- f * I_a1
    delta_I_a2 <- f * I_a2
  
    delta_I_p1 <- f * I_p1
    delta_I_p2 <- f * I_p2
    
    delta_T_1 <- f * T_1
    delta_T_2 <- f * T_2
    
    I_a1 <- I_a1 - delta_I_a1
    I_a2 <- I_a2 - delta_I_a2
  
    I_p1 <- I_p1 - delta_I_p1
    I_p2 <- I_p2 - delta_I_p2
    
    T_1 <- T_1 - delta_T_1
    T_2 <- T_2 - delta_T_2
  }

  sub observation {
    symp_1 ~ poisson(Z_sk1)
    symp_2 ~ poisson(Z_sk2)
    non_symp ~ binomial(dN_tests,Z_prev)
  }
  
  sub proposal_parameter {
     beta_bar ~ truncated_gaussian(mean = beta_bar, std = 1, lower = 0, upper = 100)
     c_22 ~ truncated_gaussian(mean = c_22, std = 0.1, lower = 0, upper = 100)
     b_1 ~ truncated_gaussian(mean = b_1, std = 0.1, lower = 0, upper = 1)
     tau_11 ~ truncated_gaussian(mean = tau_11, std = 0.5, lower = 0, upper = 32)
     tau_22 ~ truncated_gaussian(mean = tau_22, std = 0.5, lower = 0, upper = 32)
     chi ~ truncated_gaussian(mean = chi, std = 0.1, lower = 0, upper = 1)
     theta_a ~ truncated_gaussian(mean = theta_a, std = 0.1, lower = 0, upper = 1)
     theta_p ~ truncated_gaussian(mean = theta_p, std = 0.1, lower = 0, upper = 1)
  }
}
