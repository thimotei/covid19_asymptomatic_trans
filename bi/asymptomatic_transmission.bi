model asymptomatic_transmission {

  const N = 3711
  const gamma_a = 1 / 7
  const gamma_prop = 1 / 2.4
  const gamma_s = 1 / 3.2
  const t_mu = 16
  const mu = 1
  const phi = 199 / 301
  const nu = 0.25

  state S
  state E
  state I_a
  state I_p
  state I_su
  state I_sk
  state C
  state R_s
  state R_n
  state Z_ns
  state Z_sk

  param beta_bar
  param b_2
  param chi
  param theta_a
  param theta_p
  param sigma1
  param sigma2

  obs symp
  obs no_symp

  sub transition (delta = 1.0) {

    inline dN_tests_fun = 290.076 // -
        2 * 187.84 * t_now +
        3 * 59.6627 * t_now ** 2 -
        4 * 10.8069 * t_now ** 3 +
        5 * 1.20729 * t_now ** 4 -
        6 * 0.0862395 * t_now ** 5 +
        7 * 0.00394737 * t_now ** 6 -
        8 * 0.000111774 * t_now ** 7 +
        9 * (1.77781e-6) * t_now ** 8 -
        10 * (1.21093e-8) * t_now ** 9

    inline beta_t = (t_now < t_mu) ? beta_bar : beta_bar * exp(-b_2 * (t_now - t_mu))
    inline mu_t = (t_now < t_mu) ? 0 : mu

    Z_ns <- 0
    Z_sk <- 0

    ode {
      dS/dt = -beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S
      dE/dt = beta_t * ((theta_a * I_a + theta_p * I_p + I_sk + I_su) / N) * S - nu * E
      dI_a/dt = chi * nu * E - gamma_a * I_a - dN_tests_fun / (S + E + I_a + I_p + C) * I_a
      dI_p/dt = (1 - chi) * nu * E - gamma_prop * I_p - dN_tests_fun / (S + E + I_a + I_p + C) * I_p
      dI_su/dt = (1 - phi) * gamma_prop * I_p - gamma_s * I_su - mu_t * I_su
      dI_sk/dt = phi * gamma_prop * I_p - gamma_s * I_sk - mu_t * I_sk; 
      dC/dt = gamma_a * I_a + gamma_s * I_sk + gamma_s * I_su
      dR_s/dt = mu_t * I_sk + mu_t * I_su
      dR_n/dt = dN_tests_fun / (S + E + I_a + I_p + C) * I_a + dN_tests_fun / (S + E + I_a + I_p + C) * I_p
      dZ_sk/dt = phi * gamma_prop * I_p
      dZ_ns/dt = dN_tests_fun / (S + E + I_a + I_p + C) * I_a + dN_tests_fun / (S + E + I_a + I_p + C) * I_p
     }
   }

  sub parameter {
    beta_bar ~ truncated_normal(2.2, 1, lower=0, upper = 5)
    b_2 ~ truncated_normal(0, 1, lower = 0)
    chi ~ uniform(0, 1)
    theta_a ~ truncated_normal(0, 1, lower = 0, upper = 1)
    theta_p ~ truncated_normal(0, 1, lower = 0, upper = 1)
    sigma1 ~ truncated_normal(0, 1, lower = 0)
    sigma2 ~ truncated_normal(0, 1, lower = 0)
  }

  sub observation {
    symp ~ normal(Z_sk, sigma1)
    no_symp ~ normal(Z_ns, sigma2)
  }
}
