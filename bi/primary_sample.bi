model primary_sample {
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
  
  const S0_1 = N_1 
  const S0_2 = N_2-1
  
  param beta_bar
  param c_22
  param b_1
  param tau_11
  param tau_22
  param chi
  param theta_a
  param theta_p
  
  input dN_tests

  state S_1      
  state S_2      
  state E_1        
  state E_2      
  state E_1e     
  state E_2e     
  state I_a1     
  state I_a2    
  state I_p1     
  state I_p2     
  state I_su1    
  state I_su2    
  state I_sk1    
  state I_sk2    
  state T_s1
  state T_s2
  state T_a1
  state T_a2
  state C_s1
  state C_s2
  state C_a1
  state C_a2
  state R_s1     
  state R_s2     
  state R_n1     
  state R_n2     
  state Z_sk1    
  state Z_sk2    
  state Z_n1     
  state Z_n2     
  
  state f           
  state delta_I_a1   
  state delta_I_a2   
  state delta_I_p1   
  state delta_I_p2   
  state delta_T_s1    
  state delta_T_s2
  state delta_T_a1    
  state delta_T_a2

  state S            
  state E            
  state Ee           
  state I_a          
  state I_p          
  state I_su         
  state I_sk         
  state T_s
  state T_a
  state C_s 
  state C_a
  state R_s          
  state R_n          

  state Z_sk   
  state Z_n   
  state Z_prev 
  state test_pos 
  
  state N_tests   
  
  state dE_a1  
  state dE_a2  
  state dE_p1  
  state dE_p2  
  state dE_s1  
  state dE_s2  
  state dE_a    
  state dE_p    
  state dE_s    
  state dE_tot  
  state E_a    
  state E_p    
  state E_s    
  state E_tot  
  state P_a  
  state P_p  
  state P_s  

  state I_total 
  state R_total 
  state per_inf
  state per_found
  state per_missed
  state per_exp 
  state per_not_test
  state per_rec 
  state num_p_and_a
  state num_s
  state per_found_p_and_a
  state per_found_s

  state f_15 
  state f_16 
  state f_17 
  state f_18 
  state f_19 
  state f_110 
  state f_111 
  state f_112 
  state f_25 
  state f_26 
  state f_27 
  state f_28 
  state f_29 
  state f_210
  state f_211
  state f_212 
    
  state v_11
  state v_22
  state v_31
  state v_33
  state v_42
  state v_44
  state v_53
  state v_55
  state v_64
  state v_66
  state v_73
  state v_77
  state v_84
  state v_88
  state v_97
  state v_99 
  state v_108 
  state v_1010
  state v_117 
  state v_1111
  state v_128 
  state v_1212
  
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
    T_s1 <- 0
    T_s2 <- 0
    T_a1 <- 0
    T_a2 <- 0  
    C_s1 <- 0
    C_s2 <- 0
    C_a1 <- 0
    C_a2 <- 0
    R_s1 <- 0
    R_s2 <- 0
    R_n1 <- 0
    R_n2 <- 0
    
    Z_sk1 <- 0
    Z_sk2 <- 0
    Z_n1 <- 0
    Z_n2 <- 0
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
    
    Z_n1 <- (t_now % 1 == 0 ? 0 : Z_n1)
    Z_n2 <- (t_now % 1 == 0 ? 0 : Z_n2)
  
    dE_a1 <- (t_now % 1 == 0 ? 0 : dE_a1)
    dE_a2 <- (t_now % 1 == 0 ? 0 : dE_a2)
    dE_p1 <- (t_now % 1 == 0 ? 0 : dE_p1)
    dE_p2 <- (t_now % 1 == 0 ? 0 : dE_p2)
    dE_s1 <- (t_now % 1 == 0 ? 0 : dE_s1)
    dE_s2 <- (t_now % 1 == 0 ? 0 : dE_s2)
    dE_a <- (t_now % 1 == 0 ? 0 : dE_a)
    dE_p <- (t_now % 1 == 0 ? 0 : dE_p)
    dE_s <- (t_now % 1 == 0 ? 0 : dE_s)
    dE_tot <- (t_now % 1 == 0 ? 0 : dE_tot)
    
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
    
      dT_s1/dt = alpha * gamma_s_t * I_sk1 + alpha * gamma_s_t * I_su1 - eta * T_s1
      dT_s2/dt = alpha * gamma_s_t * I_sk2 + alpha * gamma_s_t * I_su2 - eta * T_s2
      
      dT_a1/dt = alpha * gamma_a * I_a1 - eta * T_a1
      dT_a2/dt = alpha * gamma_a * I_a2 - eta * T_a2
      
      dC_s1/dt = eta * T_s1 + (1 - alpha) * gamma_s_t * I_sk1 + (1 - alpha) * gamma_s_t * I_su1
      dC_s2/dt = eta * T_s2 + (1 - alpha) * gamma_s_t * I_sk2 + (1 - alpha) * gamma_s_t * I_su2
      
      dC_a1/dt = (1 - alpha) * gamma_a * I_a1 + eta * T_a1  
      dC_a2/dt = (1 - alpha) * gamma_a * I_a2 + eta * T_a2 
      
      dR_s1/dt = mu_t * I_sk1 + mu_t * I_su1
      dR_s2/dt = mu_t * I_sk2 + mu_t * I_su2
      
      dR_n1/dt = 0
      dR_n2/dt = 0
      
      dZ_sk1/dt = phi * gamma_p * I_p1
      dZ_sk2/dt = phi * gamma_p * I_p2
    
      dZ_n1/dt = 0
      dZ_n2/dt = 0
      
      ddE_a1/dt = beta_t11 * c_11 * ((theta_a * I_a1) / N_1) * S_1 + beta_t12 * c_12 * ((theta_a * I_a2) / N_2) * S_1 
      ddE_a2/dt = beta_t12 * c_12 * ((theta_a * I_a1) / N_1) * S_2 + beta_t22 * c_22 * ((theta_a * I_a2) / N_2) * S_2 
      
      ddE_p1/dt = beta_t11 * c_11 * ((theta_p * I_p1) / N_1) * S_1 + beta_t12 * c_12 * ((theta_p * I_p2) / N_2) * S_1 
      ddE_p2/dt = beta_t12 * c_12 * ((theta_p * I_p1) / N_1) * S_2 + beta_t22 * c_22 * ((theta_p * I_p2) / N_2) * S_2
      
      ddE_s1/dt = beta_t11 * c_11 * ((I_sk1 + I_su1) / N_1) * S_1 + beta_t12 * c_12 * ((I_sk2 + I_su2) / N_2) * S_1
      ddE_s2/dt = beta_t12 * c_12 * ((I_sk1 + I_su1) / N_1) * S_2 + beta_t22 * c_22 * ((I_sk2 + I_su2) / N_2) * S_2 
    }
    
    dE_a <- dE_a1 + dE_a2
    dE_p <- dE_p1 + dE_p2
    dE_s <- dE_s1 + dE_s2
    dE_tot <- dE_a + dE_p + dE_s
      
    E_a <- E_a + dE_a
    E_p <- E_p + dE_p
    E_s <- E_s + dE_s
    E_tot <- E_tot + dE_tot 
    
    P_a <- E_a/E_tot
    P_p <- E_p/E_tot
    P_s <- E_s/E_tot
    
    Z_prev <- (I_a1 + I_a2 + I_p1 + I_p2 + T_s1 + T_s2 + T_a1 + T_a2) / (S_1 + S_2 + E_1 + E_2 + E_1e + E_2e + I_a1 + I_a2 + I_p1 + I_p2 + T_s1 + T_s2 + T_a1 + T_a2 + C_s1 + C_s2 + C_a1 + + C_a2)
    
    test_pos <- dN_tests * Z_prev
    
    f <- dN_tests / (S_1 + S_2 + E_1 + E_2 + E_1e + E_2e + I_a1 + I_a2 + I_p1 + I_p2 + T_s1 + T_s2 + T_a1 + T_a2 + C_s1 + C_s2 + C_a1 + C_a2)
  
    delta_I_a1 <- f * I_a1
    delta_I_a2 <- f * I_a2
  
    delta_I_p1 <- f * I_p1
    delta_I_p2 <- f * I_p2
    
    delta_T_s1 <- f * T_s1
    delta_T_s2 <- f * T_s2
    
    delta_T_a1 <- f * T_a1
    delta_T_a2 <- f * T_a2
    
    I_a1 <- I_a1 - delta_I_a1
    I_a2 <- I_a2 - delta_I_a2
  
    I_p1 <- I_p1 - delta_I_p1
    I_p2 <- I_p2 - delta_I_p2
    
    T_s1 <- T_s1 - delta_T_s1
    T_s2 <- T_s2 - delta_T_s2
    
    T_a1 <- T_a1 - delta_T_a1
    T_a2 <- T_a2 - delta_T_a2
  
    R_n1 <- R_n1 + delta_I_a1 + delta_I_p1 + delta_T_s1 + delta_T_a1
    R_n2 <- R_n2 + delta_I_a2 + delta_I_p2 + delta_T_s2 + delta_T_a2
  
    Z_n1 <- Z_n1 + delta_I_a1 + delta_I_p1 + delta_T_s1 + delta_T_a1
    Z_n2 <- Z_n2 + delta_I_a2 + delta_I_p2 + delta_T_s2 + delta_T_a2
    
    S <- S_1 + S_2
    E <- E_1 + E_2
    Ee <- E_1e + E_2e
    I_a <- I_a1 + I_a2
    I_p <- I_p1 + I_p2
    I_su <- I_su1 + I_su2
    I_sk <- I_sk1 + I_sk2
    T_s <- T_s1 + T_s2
    T_a <- T_a1 + T_a2
    C_s <- C_s1 + C_s2
    C_a <- C_a1 + C_a2
    R_s <- R_s1 + R_s2
    R_n <- R_n1 + R_n2
    Z_sk <- Z_sk1 + Z_sk2
    Z_n <- Z_n1 + Z_n2
    
    I_total <- N_1 + N_2 - S
    R_total <- R_n + R_s
    
    per_inf <- I_total/(N_1 + N_2) 
    per_found <- R_total/I_total
    per_missed <- (I_total-R_total)/I_total
    per_exp <- (E + Ee)/I_total
    per_not_test <- (I_a + I_p + T_a + T_s + I_sk + I_su)/I_total
    per_rec <- (C_a + C_s)/I_total

    num_p_and_a <- I_p + I_a + T_a + C_a + R_n
    num_s <- I_sk + I_su + T_s + C_s + R_s
    per_found_p_and_a <- R_n/num_p_and_a
    per_found_s <- R_s/num_s
  
    N_tests <- N_tests + dN_tests 
    
    f_15 <- beta_t11 * c_11 * theta_a * S0_1 / N_1
    f_16 <- beta_t12 * c_12 * theta_a * S0_1 / N_2
    f_25 <- beta_t12 * c_12 * theta_a * S0_2 / N_1
    f_26 <- beta_t22 * c_22 * theta_a * S0_2 / N_2
    
    f_17 <- beta_t11 * c_11 * theta_p * S0_1 / N_1
    f_18 <- beta_t12 * c_12 * theta_p * S0_1 / N_2
    f_27 <- beta_t12 * c_12 * theta_p * S0_2 / N_1
    f_28 <- beta_t22 * c_22 * theta_p * S0_2 / N_2
    
    f_19 <- beta_t11 * c_11 * S0_1 / N_1
   f_110 <- beta_t12 * c_12 * S0_1 / N_2
    f_29 <- beta_t12 * c_12 * S0_2 / N_1
   f_210 <- beta_t22 * c_22 * S0_2 / N_2
    
   f_111 <- beta_t11 * c_11 * S0_1 / N_1
   f_112 <- beta_t12 * c_12 * S0_1 / N_2
   f_211 <- beta_t12 * c_12 * S0_2 / N_1
   f_212 <- beta_t22 * c_22 * S0_2 / N_2
    
      v_11 <- nu
      v_22 <- nu
      v_31 <- - nu
      v_42 <- - nu
      v_33 <- nu                   
      v_44 <- nu                   
      v_73 <- -(1-chi) * nu
      v_84 <- -(1-chi) * nu
      v_77 <- gamma_p              
      v_88 <- gamma_p              
      v_97 <- - phi * gamma_p
     v_108 <- - phi * gamma_p
     v_117 <- - (1-phi) * gamma_p
     v_128 <- - (1-phi) * gamma_p
      v_99 <- gamma_s_t + mu_t
    v_1010 <- gamma_s_t + mu_t
    v_1111 <- gamma_s_t + mu_t
    v_1212 <- gamma_s_t + mu_t
      v_53 <- - chi * nu  
      v_64 <- - chi * nu
      v_55 <- gamma_a   
      v_66 <- gamma_a   
      
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
     tau_11 ~ truncated_gaussian(mean = tau_11, std = 1, lower = 0, upper = 32)
     tau_22 ~ truncated_gaussian(mean = tau_22, std = 1, lower = 0, upper = 32)
     chi ~ truncated_gaussian(mean = chi, std = 0.1, lower = 0, upper = 1)
     theta_a ~ truncated_gaussian(mean = theta_a, std = 0.1, lower = 0, upper = 1)
     theta_p ~ truncated_gaussian(mean = theta_p, std = 0.1, lower = 0, upper = 1)
  }
}
