library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")

### Inputs ###

# Scenario name for saving the posterior and results 

posterior_name = "primary"

# LibBi model file to run MCMC on 

libbi_file_name = "primary"

# Number of iterations in MCMC 

n_samples = 1500000

### Data ###

# Symptomatic - Crew

data_cases_symp_1 <- 
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk1 = crew) %>%
  dplyr::select(time = T, value = I_sk1)

# Symptomatic - Passengers

data_cases_symp_2 <- 
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk2 = pass) %>%
  dplyr::select(time = T, value = I_sk2)

# Non-symptomatic - All 

data_cases_non_symp <- 
  readr::read_csv(file = here::here("data", "data_cases_non_symp.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>% 
  dplyr::select(time = T, value = R_n)

# Testing

data_tests <- 
  readr::read_csv(file = here::here("data", "data_tests.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, dN_tests = tests_non_symp) %>% 
  dplyr::select(time = T, value = dN_tests)

# Format data for LibBi 

obs <- list(symp_1 = data_cases_symp_1, symp_2 = data_cases_symp_2, non_symp = data_cases_non_symp)  

input <- list(dN_tests = data_tests)

### Sampling ###

model <- rbi::bi_model(file = here::here("bi", paste0(libbi_file_name, ".bi")))

# Trial fit

initial_fit = rbi::sample(model, target = "posterior", nsamples = 1000, nparticles = 1,
                          obs = obs, input = input, end_time = 32, noutputs = 32,
                          verbose = TRUE)

# Adapt proposal distributions

adapted <- initial_fit %>%
  adapt_proposal(min = 0.2, max = 0.3, adapt = "both", max_iter = 25, truncate = TRUE)

# Run MCMC

posterior <- adapted %>%
  sample(nsamples = n_samples, verbose= TRUE)

### Save posterior ###

save_libbi(posterior, here::here("posteriors", paste0(posterior_name,".rds")))


