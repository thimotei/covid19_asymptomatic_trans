library("patchwork")

here::here() %>% setwd()

source("covid19_asymptomatic_trans/R/plotting_functions.R")

posterior <- rbi::read_libbi("covid19_asymptomatic_trans/results/posterior.rds")
trajectories = summary(posterior, type = "state", quantiles = c(0.025,0.975))


breaksPlot <- seq(0, 32, 9)
datesPlot <- c("20th Jan", "29th Jan", "7th Feb", "16th Feb")
non_obs <- c("S", "E", "I_p", "I_sk", "I_su", "I_a", "C", "R_s", "R_n", "N_tests") 
titles <- c("S", "E", expression(I[p]),  expression(I[sk]),  expression(I[su]),  expression(I[a]), "C",  expression(R[s]),  expression(R[n]),  expression(N[tests])) 

plots <- list()
for(i in 1:length(non_obs_plot))
{
  colors <- seq(1, 200, length(non_obs_plot))
  p <- plottingFun(non_obs_plot[i], colors[i], titles[i], 50, withData = FALSE)
  plots[[i]] = p
}

(plots[[1]] | plots[[3]] | plots[[5]] | plots[[7]] | plots[[9]])/
(plots[[2]] | plots[[4]] | plots[[6]] | plots[[8]] | plots[[10]])
