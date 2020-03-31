parametersToPlot <- c("beta_bar", "b_2", "chi", "theta_p", "sigma1", "sigma2")
draws <- as.array(fit_full, pars=parametersToPlot)
bayesplot::mcmc_trace(draws)

posterior <- as.matrix(fit_full)
plot_title <- ggplot2::ggtitle("Posterior distributions",
                               "with medians and 80% intervals")
bayesplot::mcmc_areas(posterior,
                      pars = parametersToPlot) + plot_title

pairs(fit_full, pars = c("beta_bar", "sigma1"))

pairs(fit_full, pars = parametersToPlot)