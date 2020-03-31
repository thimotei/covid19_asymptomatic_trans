parametersToPlot <- c("beta_bar", "b_2", "chi", "theta_a")
draws <- as.array(fit_full, pars=parametersToPlot)
bayesplot::mcmc_trace(draws)

posterior <- as.matrix(fit_full)
plot_title <- ggplot2::ggtitle("Posterior distributions",
                               "with medians and 80% intervals")
bayesplot::mcmc_areas(posterior,
                      pars = parametersToPlot) + plot_title
