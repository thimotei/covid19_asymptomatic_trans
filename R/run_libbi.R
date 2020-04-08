#library('devtools')
#install_github("libbi/rbihelpers")
#library(devtools)
#install_github("sbfnk/rbi.helpers")

library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")
library("GGally")
library('coda')

here::here()

data_cases_symp <- 
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk = all) %>%
  dplyr::select(time = T, value = I_sk)

data_cases_non_symp <- 
  readr::read_csv(file = here::here("data", "data_cases_non_symp.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>% 
  dplyr::select(time = T, value = R_n)

data_tests <- 
  readr::read_csv(file = here::here("data", "data_tests.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, dN_tests = tests_non_symp) %>% 
  dplyr::select(time = T, value = dN_tests)

obs <- list(symp = data_cases_symp, non_symp = data_cases_non_symp)

input <- list(dN_tests = data_tests)

# init <- list(beta_bar = 2.0,
#              b_1 = 0.95,
#              b_2 = 10,
#              tau = 19,
#              chi = 0.825,
#              theta_a = 0.25,
#              theta_p = 0.5)

model <- rbi::bi_model(file = here::here("bi", "asymptomatic_transmission_post_infection.bi"))               
model <- fix(model, b_2 = 10, alpha = 0)                                                            
#model <- fix(model, b_2 = 10, chi = 0, theta_a = 0)  
 
initial_fit = rbi::sample(model, target = "posterior", nsamples = 10000, nparticles = 1,
                          obs = obs, input = input, end_time = 32, noutputs = 32,
                          verbose = TRUE)

adapted <- initial_fit %>%
  adapt_proposal(min = 0.2, max = 0.3)
 

posterior <- adapted %>%
  sample(nsamples = 500000)

#save_libbi(posterior, here::here("results", "posterior_baseline_long.rds"))

#bi_read(posterior,thin = 10)






