// ode_asymptomatic transmission model
functions {
  real beta_t(real beta_bar, real b_1, real b_2, real tau, real t)
  {
    return(beta_bar*(1-b_1/(1+exp(-b_2*(t-tau)))));
  }
  real mu_t(real t, real t_mu)
  {
    real mu_t;
    real mu;
    if(t < t_mu)
    {
        mu_t = 0;
    }
    else
    {
        mu_t = mu;
    }
      return(mu_t);
  }
  real dN_tests_fun(real t)
  {
    return(-1.67487+13.8278*t-6.26078*t^2+0.841799*t^3-0.0438224*t^4+0.000838706*t^5);
  }
  real[] ode_model(real t, 
                   real[] y, 
                   real[] theta, 
                   real[] y_r,
                   int[] y_i) 
  {
    
    real dydt[9];
    dydt[1] = -beta_t(theta[1], theta[2], theta[3], theta[4], t)*((theta[6]*y[3] + theta[13]*y[4] + y[8] + y[5])/theta[15])*y[1];
    dydt[2] = -beta_t(theta[1], theta[2], theta[3], theta[4], t)*((theta[6]*y[3] + theta[13]*y[4] + y[8] + y[5])/theta[15])*y[1] - theta[7]*y[2];
    dydt[3] = theta[5]*theta[7]*y[2] - theta[8]*y[3] - dN_tests_fun(t)/(y[1]+y[2]+y[4]+y[3]+y[6])*y[3];
    dydt[4] = (1-theta[5])*theta[7]*y[2] - theta[9]*y[4] - dN_tests_fun(t)/(y[1]+y[2]+y[4]+y[3]+y[6])*y[4];
    dydt[5] = (1-theta[14])*theta[9]*y[4] - theta[10]*y[5] - mu_t(t, theta[11])*y[5];
    dydt[6] = theta[8]*y[3] + theta[10]*y[8] + theta[10]*y[5];
    dydt[7] = mu_t(t, theta[11])*y[8] + mu_t(t, theta[11])*y[5];
    dydt[8] = theta[14]*theta[9]*y[4] - theta[10]*y[8] - mu_t(t, theta[11])*y[8];
    dydt[9] = dN_tests_fun(t)/(y[1]+y[2]+y[4]+y[3]+y[6])*y[3] + dN_tests_fun(t)/(y[1]+y[2]+y[4]+y[3]+y[6])*y[4];
    return dydt;
  }
}
data {
  int<lower=1> t;
  real<lower=0> y_obs_1[t];
  real<lower=0> y_obs_2[t];
  real t0;
  real y0[9];
  real ts[t];
}
transformed data {
  real y_r[0];
  int y_i[0];
}
parameters {
  real<lower=0> sigmaNoise1;
  real<lower=0> sigmaNoise2;
  real<lower=0> beta_bar;
  real<lower=0> b_1;         
  real<lower=0> b_2;        
  real<lower=0> tau;    
  real<lower=0> chi; 
  real<lower=0> theta_a;    
  real<lower=0> y_mis[t,7];
}
transformed parameters{
    real y_hat[t,9];
    {
      
      real nu = 1.0/4.0;          
      real gamma_a = 1.0/7.0;     
      real gamma_prop = 1.0/2.4;   
      real gamma_s = 1.0/3.2;   
      real t_mu = 16.0;         
      real mu = 1.0;            
      real theta_p = 0.99;
      real phi = 199.0/301.0;
      int N = 3711;
      
      real theta[15];
      theta[1] = beta_bar;
      theta[2] = b_1;
      theta[3] = b_2;
      theta[4] = tau;
      theta[5] = chi;
      theta[6] = theta_a;
      theta[7] = nu;
      theta[8] = gamma_a;
      theta[9] = gamma_prop;
      theta[10] = gamma_s;
      theta[11] = t_mu;
      theta[12] = mu;
      theta[13] = theta_p;
      theta[14] = phi;
      theta[15] = N;
    
      y_hat = integrate_ode_bdf(ode_model, y0, t0, ts, theta, y_r, y_i);

    }
}

model {

  beta_bar ~ uniform(0,2);
  b_1 ~ uniform(0,1);         
  b_2 ~ uniform(0,1000);        
  tau ~ uniform(0,32);
  chi ~ uniform(0,1);
  theta_a ~ uniform(0,1);
  
  sigmaNoise1 ~ cauchy(0,2.5);
  sigmaNoise2 ~ cauchy(0,2.5);
  
  for (i in 1:t)
  {
    y_mis[i,7] ~ normal(y_hat[i, 1:7], 0.01);
    y_obs_1[i] ~ normal(y_hat[i, 8], sigmaNoise1);
    y_obs_2[i] ~ normal(y_hat[i, 9], sigmaNoise2);
    
  }
}
    