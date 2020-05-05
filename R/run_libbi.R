library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")

### Inputs ###

posterior_name = ""

libbi_file_name = ""

model <- rbi::bi_model(file = here::here("bi", paste0(libbi_file_name, ".bi")))
model <- fix(model)

n_samples = 2000000

### Data ###

data_cases_symp_1 <- 
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk1 = crew) %>%
  dplyr::select(time = T, value = I_sk1)

data_cases_symp_2 <- 
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk2 = pass) %>%
  dplyr::select(time = T, value = I_sk2)

data_cases_non_symp <- 
  readr::read_csv(file = here::here("data", "data_cases_non_symp.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>% 
  dplyr::select(time = T, value = R_n)

data_tests <- 
  readr::read_csv(file = here::here("data", "data_tests.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, dN_tests = tests_non_symp) %>% 
  dplyr::select(time = T, value = dN_tests)

input <- list(dN_tests = data_tests)

### Sampling

initial_fit = rbi::sample(model, target = "posterior", nsamples = 1000, nparticles = 1,
                          obs = obs, input = input, end_time = 32, noutputs = 32, init = init,
                          verbose = TRUE)

adapted <- initial_fit %>%
  adapt_proposal(min = 0.2, max = 0.3, adapt = "both", max_iter = 25, truncate = TRUE)

posterior <- adapted %>%
  sample(nsamples = n_samples, verbose= TRUE)

save_libbi(posterior, here::here("posteriors", paste0(posterior_name,".rds")))


