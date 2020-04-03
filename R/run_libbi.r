library("tidyverse")
library("rbi")
library("rbi.helpers")
library("GGally")

here::here()

data_cases_symp <-
  readr::read_csv(file = here::here("data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk = all) %>%
  dplyr::select(time = T, value = I_sk)

data_cases_non_symp <-
  readr::read_csv(file = here::here("data", "data_cases_non_symp.csv") %>%
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>%
  dplyr::select(time = T, value = R_n)

obs <- list(symp = data_cases_symp, non_symp = data_cases_non_symp)

init <- list(S = 3710.0,
             E = 0,
             I_A = 0,
             I_p = 0,
             I_su = 0,
             I_sk = 1,
             C = 0,
             R_s = 0,
             R_n = 0,
             Z_sk = 0,
             Z_ns = 0)

model <-
  rbi::bi_model(file = here::here("bi", "asymptomatic_transmission.bi")

initial_fit <-
  rbi::sample(model, target = "posterior", proposal = "prior", nsamples = 1000,
              init = init, obs = obs, end_time = 32, noutputs = 32,
              verbose = TRUE)

adapted <- initial_fit %>%
  adapt_proposal(min = 0.2, max = 0.3)

save_libbi(adapted, here::here("results", "adapted.rds"))

posterior <- adapted %>%
  sample(nsamples = 10000)

traces <- get_traces(posterior, thin = 10)

## pairs plot
ggpairs(traces)

