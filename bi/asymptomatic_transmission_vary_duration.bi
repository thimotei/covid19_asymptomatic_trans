model asymptomatic_transmission {
  const N = 3711
  const nu = 1 / 4
  //const gamma_a = 1 / 7
  const gamma_p = 1 / 2.4
  const gamma_s = 1 / 3.2
  const eta = 1 / 7
  const t_mu = 16
  const mu = 1
  const phi = 199 / 314
  const alpha = 0
  
  param beta_bar
  param b_1
  param b_2
  param tau
  param chi
  param theta_a
  param theta_p
  param gamma_a

  input dN_tests

  state S
  state E
  state I_a
  state I_p
  state I_su
  state I_sk
  state T
  state C
  state R_s
  state R_n
  state Z_sk
  state Z_n
  state E_tot
  state E_a
  state E_p
  state E_s
  state P_a
  state P_p
  state P_s
  state dE_tot
  state dE_a
  state dE_p
  state dE_s
  state dP_a
  state dP_p
  state dP_s
  state N_tests
  state delta_I_a
  state delta_I_p
  state delta_T
  obs symp
  obs non_symp
  
  sub initial {
    S <- 3710.0
    E <- 0
    I_a <- 0
    I_p <- 0
    I_su <- 0
    I_sk <- 1
    T <- 0
    C <- 0
    R_s <- 0
    R_n <- 0
    E_tot <- 0
    E_a <- 0
    E_p <- 0
    E_s <- 0
    P_a <- 0
    P_p <- 0
    P_s <- 0
    dE_tot <- 0
    dE_a <- 0
    dE_p <- 0
    dE_s <- 0
    dP_a <- 0
    dP_p <- 0
    dP_s <- 0
    Z_sk <- 0
    Z_n <- 0
    delta_I_a <- 0
    delta_I_a <- 0  
    delta_T <- 0
    N_tests <- 0
  }

  sub parameter {
    beta_bar ~ uniform(0,100)
    b_1 ~ uniform(0, 1)
    b_2 ~ uniform(0,100)
    tau ~ uniform(0, 32)
    chi ~ uniform(0, 1)
    theta_a ~ uniform(0, 1)
    theta_p ~ uniform(0, 1)
    gamma_a ~ uniform(1/21,1)
  }

  sub transition (delta = 1.0) {

    inline beta_t = beta_bar*(1-b_1/(1+exp(-b_2*(t_now-tau))))
    inline mu_t = (t_now < t_mu) ? 0 : mu

    Z_sk <- (t_now % 1 == 0 ? 0 : Z_sk)
    Z_n <- (t_now % 1 == 0 ? 0 : Z_n)
    
    dE_tot <- (t_now % 1 == 0 ? 0 : dE_tot)
    dE_a <- (t_now % 1 == 0 ? 0 : dE_a)
    dE_p <- (t_now % 1 == 0 ? 0 : dE_p)
    dE_s <- (t_now % 1 == 0 ? 0 : dE_s)
    
    dP_a <- (t_now % 1 == 0 ? 0 : dP_a)
    dP_p <- (t_now % 1 == 0 ? 0 : dP_p)
    dP_s <- (t_now % 1 == 0 ? 0 : dP_s)
    
    ode {
      dS/dt = -beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S
      dE/dt = beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S - nu * E
      dI_a/dt = chi * nu * E - gamma_a * I_a
      dI_p/dt = (1 - chi) * nu * E - gamma_p * I_p
      dI_su/dt = (1 - phi) * gamma_p * I_p - gamma_s * I_su - mu_t * I_su
      dI_sk/dt = phi * gamma_p * I_p - gamma_s * I_sk - mu_t * I_sk
      dT/dt = alpha * gamma_a * I_a - eta * T
      dC/dt = (1 - alpha) * gamma_a * I_a + eta * T + gamma_s * I_sk + gamma_s * I_su
      dR_s/dt = mu_t * I_sk + mu_t * I_su
      dR_n/dt = 0
      dZ_sk/dt = phi * gamma_p * I_p
      dZ_n/dt = 0
      dE_tot/dt = beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S 
      dE_a/dt = beta_t * ((theta_a * I_a) / N) * S
      dE_p/dt = beta_t * ((theta_p * I_p) / N) * S
      dE_s/dt = beta_t * ((I_sk + I_su) / N) * S
      ddE_tot/dt = beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S 
      ddE_a/dt = beta_t * ((theta_a * I_a) / N) * S
      ddE_p/dt = beta_t * ((theta_p * I_p) / N) * S
      ddE_s/dt = beta_t * ((I_sk + I_su) / N) * S
    }
    
    P_a <- E_a/E_tot
    P_p <- E_p/E_tot
    P_s <- E_s/E_tot
    
    dP_a <- dE_a/dE_tot
    dP_p <- dE_p/dE_tot
    dP_s <- dE_s/dE_tot
  
    delta_I_a <- dN_tests / (S + E + I_a + I_p + T + C) * I_a
    delta_I_p <- dN_tests / (S + E + I_a + I_p + T + C) * I_p
    delta_T <- dN_tests / (S + E + I_a + I_p + T + C) * T
    
    I_a <- I_a - delta_I_a
    I_p <- I_p - delta_I_p
    T <- T - delta_T
    R_n <- R_n + delta_I_a + delta_I_p + delta_T
    Z_n <- Z_n + delta_I_a + delta_I_p + delta_T
    N_tests <- N_tests + dN_tests 
    
  }

  sub observation {
    symp ~ poisson(Z_sk)
    non_symp ~ poisson(Z_n) 
  }
  
  sub proposal_parameter {
     beta_bar ~ truncated_gaussian(mean = beta_bar, std = 1, lower = 0, upper = 100)
     b_1 ~ truncated_gaussian(mean = b_1, std = 0.1, lower = 0, upper = 1)
     b_2 ~ truncated_gaussian(mean = b_2, std = 1, lower = 0, upper = 100)
     tau ~ truncated_gaussian(mean = tau, std = 1, lower = 0, upper = 32)
     chi ~ truncated_gaussian(mean = chi, std = 0.1, lower = 0, upper = 1)
     theta_a ~ truncated_gaussian(mean = theta_a, std = 0.1, lower = 0, upper = 1)
     theta_p ~ truncated_gaussian(mean = theta_p, std = 0.1, lower = 0, upper = 1)
     gamma_a ~ truncated_gaussian(mean = gamma_a, std = 1, lower = 1/21, upper = 1)
  }
}
