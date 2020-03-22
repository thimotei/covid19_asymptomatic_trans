setwd("/Users/jonemery/Google Drive/LSHTM - Research Fellow/03 R/Working directory")
library("deSolve", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

### times ###

# Allocate all daily values (e.g. symptomatic cases per day) to the end of the day 
# E.g. 2 cases on 20th Jan (the first day) are allocated to 23:59:59 on 20th Jan (time = 1 day) 

initial_time = 0 # 00:00:00 on 20th Jan
final_time = 32  # 23:59:59 on 20th Feb 
time_step = 0.1  # days

### parameters ###

params = c(beta_bar = 1,    # 1/days      # initial transmission rate (approx) 
           b_1 = 0,       # number      # final transmission rate = (1-b_1)*initial transmission rate 
           b_2 = 5,         # number      # gradient of transition between intial and final transmission rate 
           tau = 1,         # days        # timescale of transition between intial and final transmission rate 
           nu = 1/3,        # 1/days      # 1/latent period 
           gamma_a = 1/11,  # 1/days      # 1/infectious period of asymptomatic
           gamma_p = 1/4,   # 1/days      # 1/infectious period of presymptomatic 
           gamma_s = 1/7,   # 1/days      # 1/infectious period of symptomatic 
           t_mu = 16,       # days        # time after which passive case finding starts (16 = 00:00:00 on 5th Feb)
           mu = 1,          # 1/days      # 1/delay between symptom onset and removal through passive case finding
           f = 1/3,         # 1/days      # PCR screening rate 
           theta_p = 1,     # proportion  # relative infectiousness of presymptomatic relative to symptomatic 
           theta_a = 1,     # proportion  # relative infectiousness of asymptomatic relative to symptomatic 
           chi = 0.5,       # proportion  # proportion exposed that become asymptomatic  
           phi = 199/301,   # proportion  # proportion of presymptomatic that become symptomatic with a known onset date
           N = 3711)        # number      # initial number on the ship (assumed the number as at 5th Feb applies on 20th Jan)

### model ###

traj = function(parameters){

  times = seq(initial_time,final_time,time_step)
  
  #variables and initial conditions 
  
  state = with(as.list(c(parameters)),{
    
    c(S = N-1,    # Susceptible 
      E = 0,      # Exposed 
      I_a = 0,      # Infectious (symptomatic)
      I_p = 0,      # Infectious (presymptomatic)
      I_sk = 1,      # Infectious (symptomatic - known date of onset)
      I_su = 0,      # Infectious (symptomatic - unknown date of onset)
      C = 0,      # Self-cured on the ship 
      R_s = 0,      # Symptomatics removed from the ship through passive case finding
      R_n = 0)      # Non-symptomtics (pre or asymptomatics) removed from the ship through active case finding
  })
    
  #function for the sigmoid transmission rate to reflect reduced contact after quarantine 
  
  beta_t = function(parameters,time){
    with(as.list(c(parameters)),{
      beta_t = beta_bar*(1-b_1/(1+exp(-b_2*(time-tau))))
      return(beta_t)
    })
  }
  
  #function for passive case finding of symptomatics after quarantine starts 
  
  mu_t = function(parameters,time){
    with(as.list(c(parameters)),{
      if(time<t_mu){
        mu_t = 0
      }else{
        mu_t = mu
      }
      return(mu_t)
    })
  }

  # model equations 
  
  equations = function(time, state, parameters) {
    
    with(as.list(c(state,parameters)),{
      
      dS = -beta_t(parameters,time)*((theta_a*I_a + theta_p*I_p + I_sk + I_su)/N)*S
      dE =  beta_t(parameters,time)*((theta_a*I_a + theta_p*I_p + I_sk + I_su)/N)*S - nu*E
    dI_a =  chi*nu*E - gamma_a*I_a - f*I_a
    dI_p =  (1-chi)*nu*E - gamma_p*I_p - f*I_p
   dI_sk =  phi*gamma_p*I_p - gamma_s*I_sk - mu_t(parameters,time)*I_sk
   dI_su =  (1-phi)*gamma_p*I_p - gamma_s*I_su - mu_t(parameters,time)*I_su 
      dC =  gamma_a*I_a + gamma_s*I_sk + gamma_s*I_su
    dR_s =  mu_t(parameters,time)*I_sk + mu_t(parameters,time)*I_su
    dR_n =  f*I_a + f*I_p
        
      return(list(c(dS,dE,dI_a,dI_p,dI_sk,dI_su,dC,dR_s,dR_n),
                  I_s=I_sk+I_su,                                # total symptomatics 
                  R=R_s+R_n,                                    # total removed from the ship
                  N_ship=S+E+I_a+I_p+I_sk+I_su+C,               # total remaining on the ship
                  test=S+E+I_a+I_p+I_sk+I_su+C+R_s+R_n))        # test to make sure total people remains constant 
    })
  }
  
  output = ode(func = equations, y = state, parms = parameters, times = times)
  
  return(output)
}

plot(traj(params))
