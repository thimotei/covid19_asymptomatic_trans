setwd("/Users/jonemery/Google Drive/LSHTM - Research Fellow/03 R/Working directory")
library("deSolve", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("numDeriv", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("FME", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

### times ###

# Allocate all daily values (e.g. symptomatic cases per day) to the end of the day 
# E.g. 2 cases on 20th Jan (the first day) are allocated to 23:59:59 on 20th Jan (time = 1 day) 

initial_time = 0 # 00:00:00 on 20th Jan
final_time = 32  # 23:59:59 on 20th Feb 
time_step = 0.1  # days

### parameters ###

params = c(beta_bar = 2.2,    # 1/days       - initial transmission rate (approx) 
           b_1 = 0.9,         # number       - final transmission rate = (1-b_1)*initial transmission rate 
           b_2 = 1.38,        # number       - gradient of transition between intial and final transmission rate 
           tau = 16,          # days         - timescale of transition between intial and final transmission rate 
           nu = 1/4,          # 1/days       - 1/latent period 
           gamma_a = 1/7,     # 1/days       - 1/infectious period of asymptomatic
           gamma_p = 1/2.4,   # 1/days       - 1/infectious period of presymptomatic 
           gamma_s = 1/3.2,   # 1/days       - 1/infectious period of symptomatic 
           t_mu = 16,         # days         - time after which passive case finding starts (16 = 00:00:00 on 5th Feb)
           mu = 1,            # 1/days       - 1/delay between symptom onset and removal through passive case finding
           theta_p = 0.99,    # proportion   - relative infectiousness of presymptomatic relative to symptomatic 
           theta_a = 0.5,     # proportion   - relative infectiousness of asymptomatic relative to symptomatic 
           chi = 0.3,         # proportion   - proportion exposed that become asymptomatic  
           phi = 199/301,     # proportion   - proportion of presymptomatic that become symptomatic with a known onset date
           N = 3711)          # number       - initial number on the ship (assumed the number as at 5th Feb applies on 20th Jan)

params_var = params[c("beta_bar","b_1","b_2","tau","t_mu", "theta_a", "chi", "theta_p")]  # parameters to vary in fitting
params_fix = params[c("nu","gamma_a","gamma_p","gamma_s","mu","phi","N")]                 # parameters to hold fixed in fitting

lower = log(c("beta_bar" = 0.0001, "b_1" = 0.0001, "b_2" = 0.0001, "tau" = 0.0001,  "t_mu" = 0.0001,  "theta_a" = 0.0001, "chi" = 0.0001, "theta_p" = 0.0001))  # lower bounds on parameters in fitting 
upper = log(c("beta_bar" = 1000 ,  "b_1" = 1 ,     "b_2" = 1000,   "tau" = 32,      "t_mu" = 32,      "theta_a" = 1,      "chi" = 1,      "theta_p" = 1))       # upper bounds on parameters in fitting 

### PCR screening function ###

# Screening funtion f(t) = dN_tests/dt * 1/(S+E+I_p+I_a+C). Below calculates dN_tests/dt

data_tests = read.table(file = "data_tests.csv", header = TRUE, sep = ",") # Input data on cumulative PCR tests on everyone except symptomatics and removed

N_tests = splinefun(data_tests$day_no, data_tests$tests_non_symp_cum, method = "hyman") # Convert PCR test data to a continuous function 
                                                                                        # Note that t<0 and t>32 haven't been accounted for
 
dN_tests = function(time){              # Time derivative of cumulative tests
  dN_tests = grad(N_tests, x = time)
  return(dN_tests)
}

### model ###

traj = function(parameters){

  times = seq(initial_time,final_time,time_step)
  
  #variables and initial conditions 
  
  state = with(as.list(c(parameters)),{
    
    c(S = N-1,      # Susceptible 
      E = 0,        # Exposed 
      I_a = 0,      # Infectious (symptomatic)
      I_p = 0,      # Infectious (presymptomatic)
      I_sk = 1,     # Infectious (symptomatic - known date of onset)
      I_su = 0,     # Infectious (symptomatic - unknown date of onset)
      C = 0,        # Self-cured on the ship 
      R_s = 0,      # Symptomatics removed from the ship through passive case finding
      R_n = 0,      # Non-symptomtics (pre or asymptomatics) removed from the ship through active case finding
      
      S_tests = 0,  # Number of tests performed on Susceptible 
      E_tests = 0,  # Number of tests performed on Exposed
      C_tests = 0)  # Number of tests performed on Self-cured
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
     dI_a =  chi*nu*E - gamma_a*I_a - dN_tests(time)/(S+E+I_p+I_a+C)*I_a
     dI_p =  (1-chi)*nu*E - gamma_p*I_p - dN_tests(time)/(S+E+I_p+I_a+C)*I_p
    dI_sk =  phi*gamma_p*I_p - gamma_s*I_sk - mu_t(parameters,time)*I_sk
    dI_su =  (1-phi)*gamma_p*I_p - gamma_s*I_su - mu_t(parameters,time)*I_su 
       dC =  gamma_a*I_a + gamma_s*I_sk + gamma_s*I_su
     dR_s =  mu_t(parameters,time)*I_sk + mu_t(parameters,time)*I_su
     dR_n =  dN_tests(time)/(S+E+I_p+I_a+C)*I_a + dN_tests(time)/(S+E+I_p+I_a+C)*I_p
     
 dS_tests = dN_tests(time)/(S+E+I_p+I_a+C)*S  
 dE_tests = dN_tests(time)/(S+E+I_p+I_a+C)*E   
 dC_tests = dN_tests(time)/(S+E+I_p+I_a+C)*C   
        
      return(list(c(dS,dE,dI_a,dI_p,dI_sk,dI_su,dC,dR_s,dR_n,dS_tests,dE_tests,dC_tests),
                  I_s= I_sk+I_su,                                # total symptomatics 
                  R= R_s+R_n,                                    # total removed from the ship
                  N_ship= S+E+I_a+I_p+I_sk+I_su+C,               # total remaining on the ship
                  N_tests = S_tests+E_tests+C_tests+R_n,         # total tested through PCR screening 
                  test=S+E+I_a+I_p+I_sk+I_su+C+R_s+R_n,          # test to make sure total people remains constant 
                  dI_sk =  phi*gamma_p*I_p,                      # daily symptomatic cases with an onset date - for fitting 
                  dR_n = dN_tests(time)/(S+E+I_p+I_a+C)*I_a + dN_tests(time)/(S+E+I_p+I_a+C)*I_p)) # daily non-symptomatic cases removed from the ship following testing - for fitting 
    })
  }
  
  output = ode(func = equations, y = state, parms = parameters, times = times)
  
  return(output)
}

traj_fit = function(parameters){ 
  output = traj(parameters)
  print(parameters)
  return(as.data.frame(output))
}

### fitting data ###

data_cases_symp = read.csv(file = "data_cases_symp.csv", header = TRUE, sep = ",")          # bring in symptomatic cases data for fitting 
names(data_cases_symp)[c(2,3)] = c("time","dI_sk")                                          # convert names to match model outputs 
data_cases_symp$av = mean(data_cases_symp$dI_sk)                                            # add a mean column for fitting 

data_cases_non_symp = read.csv(file = "data_cases_non_symp.csv", header = TRUE, sep = ",")  # bring in non-symptomatic cases data for fitting 
names(data_cases_non_symp)[c(2,3)] = c("time","dR_n")                                       # convert names to match model outputs 
data_cases_non_symp$av = mean(data_cases_non_symp$dR_n)                                     # add a mean column for fitting 

### least squares fitting ###

# define a cost function that compares model outputs with data 

cost_function = function(parameters){
  output = traj_fit(parameters)
  cost_symp = modCost(model = output, obs = data_cases_symp[,c("time","dI_sk","av")], err = "av")
  cost_non_symp = modCost(model = output, obs = data_cases_non_symp[,c("time","dR_n","av")], cost = cost_symp, err = "av")
}

# set up cost function to input into optimisation 

log_cost_function = function(log_parameters){
  cost_function(c(exp(log_parameters), params_fix))
}

# optimisation 

params_ls = modFit(f = log_cost_function, p = log(params_var), method = "Pseudo", lower = lower, upper = upper) # log-least squares parameters  
 
params_ls$par = exp(coef(params_ls)) # unlogged least squares parameters  
 
params_ls = c(params_ls$par, params_fix) # stitch fitted parameters back togeter with fixed parameters 

# plot results 

plot(traj(params_ls))

plot(traj(params_ls)[,"time"],traj(params_ls)[,"dR_n"],type="l")
points(data_cases_non_symp$time,data_cases_non_symp$dR_n)

plot(traj(params_ls)[,"time"],traj(params_ls)[,"dI_sk"],type="l")
points(data_cases_symp$time,data_cases_symp$dI_sk)















