#setwd(here::here())

data_cases_symp <- readr::read_csv(file = "data/data_cases_symp.csv") %>% 
  dplyr::mutate(date = onset_date, T = day_no, I_sk = all) %>%
  dplyr::select(T, I_sk)

data_cases_non_symp <- readr::read_csv(file = "data/data_cases_non_symp.csv") %>% 
  dplyr::mutate(date = test_date, T = day_no, R_n = non_symp) %>%
  dplyr::select(T, R_n)

#data_together <- dplyr::left_join(data_cases_symp, data_cases_non_symp) 

data_together <- list(t = length(data_cases_symp$I_sk),
                      y_obs_1 = data_cases_symp$I_sk,
                      y_obs_2 = data_cases_non_symp$R_n,
                      t0 = 0,
                      y0 = c(3710.0,0,0,0,0,0,0,1,0),
                      ts = seq(1,32,1))


mod_full <- rstan::stan_model(file = "stan/asymptomatic_transmission.stan")

fit_full <- rstan::sampling(mod_full,
                            data = data_together, 
                            chains = 1, 
                            verbose = TRUE,
                            iter = 100,
                            cores = 4,
                            refresh = 20)
                            #control = list(adapt_delta = 0.8))

res <- rstan::extract(fit_full)
res$y_hat[,,1] # this is yhat[,1]

data.frame(t = 1:32, 
           S = res$y_hat[,,1],
           E = res$y_hat[,,2],
           I_a = res$y_hat[,,3],
           I_p = res$y_hat[,,4],
           I_su = res$y_hat[,,5],
           C = res$y_hat[,,6],
           R_s = res$y_hat[,,7], 
           I_sk = res$y_hat[,,8],
           R_n = res$y_hat[,,9]) %>% 
  tidyr::gather("compartment","value", -t) %>%
  ggplot2::ggplot(ggplot2::aes(x = t, y = value, col = compartment)) +
  ggplot2::geom_line()