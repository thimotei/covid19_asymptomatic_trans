library("patchwork")

here::here() %>% setwd()

source("covid19_asymptomatic_trans/R/plotting_functions.R")


posterior <- rbi::read_libbi("covid19_asymptomatic_trans/results/posterior.rds")
trajectories = summary(posterior, type = "state", quantiles = c(0.025,0.975))

breaksPlot <- seq(0, 32, 9)
datesPlot <- c("20th Jan", "29th Jan", "7th Feb", "16th Feb")
obs <- c("Z_n",  "Z_sk") 
titles <- c("Clinical cases", "Subclinical cases")

data_cases_symp <- 
  readr::read_csv(file = here::here("covid19_asymptomatic_trans/data", "data_cases_symp.csv")) %>%
  dplyr::mutate(date = onset_date, T = day_no, I_sk = all) %>%
  dplyr::select(time = T, value = I_sk)

data_cases_non_symp <- 
  readr::read_csv(file = here::here("covid19_asymptomatic_trans/data", "data_cases_non_symp.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>% 
  dplyr::select(time = T, value = R_n)

data_tests <- 
  readr::read_csv(file = here::here("covid19_asymptomatic_trans/data", "data_tests.csv")) %>% 
  dplyr::mutate(date = test_date, T = day_no, dN_tests = tests_non_symp) %>% 
  dplyr::select(time = T, value = dN_tests)

dates <- 
  readr::read_csv(file = here::here("covid19_asymptomatic_trans/data", "dates_day_no_conv.csv")) %>%
  dplyr::select(dates = date)

dataTogether <- data_cases_non_symp %>%
  dplyr::left_join(data_cases_symp, by = "time") %>%
  dplyr::rename(non_symp = value.x, 
                symp = value.y)

plots <- list()
for(i in 1:length(obs))
{
  colors <- seq(1, 200, length(obs))
  p <- plottingFun(obs[i], colors[i], titles[i], 100, withData = TRUE)
  plots[[i]] = p
}


plots[[2]] | plots[[1]]
