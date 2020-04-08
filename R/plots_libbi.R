library('ggplot2')
library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")
library("GGally")
library('coda')
library('bayesplot')
library("here")

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

#posterior = read_libbi(here::here("results", paste0("posterior_",scenario,".rds")))

#scenario = "baseline_long"

traces = get_traces(posterior, thin = 10)

trajectories_resample = predict(posterior, start_time=0, end_time=32, output_every=1, nsamples = 10000)

trajectories = summary(trajectories_resample, type = "state")
trajectories_non_obs = subset(trajectories, var %in% c("S","E","I_p","I_sk","I_su","I_a","T","C","R_s","R_n","N_tests"))
trajectories_obs = subset(trajectories, var %in% c("Z_sk","Z_n"))
trajectories_prop_trans = subset(trajectories, var %in% c("dP_a","dP_p","dP_s"))
trajectories_prop_trans_cum = subset(trajectories, var %in% c("P_a","P_p","P_s"))

title_size = 20
axis_title_size = 18
axis_text_size = 14
title_align = 0.5
title_face = "bold"
line_thickness = 0.7
colour_line = "#CE2931"
point_size = 2
alpha = 0.2
height = 16
aspect_ratio = 1.4

a = mcmc_pairs(traces,diag_fun = c("dens"), off_diag_args = list(size = 0.5, alpha = 0.2))

ggsave(a, file = here::here("plots", paste0(scenario,"_pairs.png")), width = 1.3*aspect_ratio*height, height = height, units = "cm", dpi = 300)

b = mcmc_trace(traces)

ggsave(b, file = here::here("plots", paste0(scenario,"_traces.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

c = mcmc_acf(traces)

ggsave(c, file = here::here("plots", paste0(scenario,"_acf.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

d = ggplot(trajectories_non_obs, aes(x=time)) +
  geom_line(aes(y=Median), colour = colour_line) +
  geom_ribbon(aes(ymin= Min., ymax= Max.), fill = colour_line, alpha=alpha) +
  facet_wrap(vars(var), scales="free",nrow=3) +
  ylab("Number") 

ggsave(d, file = here::here("plots", paste0(scenario,"_trajs.png")), width = 1.3*aspect_ratio*height, height = height, units = "cm", dpi = 300)

e = ggplot(trajectories_prop_trans, aes(x=time)) +
  geom_line(aes(y=Median), colour = colour_line) +
  geom_ribbon(aes(ymin= Min., ymax= Max.), fill = colour_line, alpha=alpha) +
  facet_wrap(vars(var), scales="fix",nrow=2) +
  ylab("Proportion") 

ggsave(e, file = here::here("plots", paste0(scenario,"_prop_trans.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

f = ggplot(trajectories_prop_trans_cum, aes(x=time)) +
  geom_line(aes(y=Median), colour = colour_line) +
  geom_ribbon(aes(ymin= Min., ymax= Max.), fill = colour_line, alpha=alpha) +
  facet_wrap(vars(var), scales="fix",nrow=2) +
  ylab("Proportion") 

ggsave(f, file = here::here("plots", paste0(scenario,"_prop_trans_cum.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

g = ggplot(subset(trajectories_obs, var == "Z_sk"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin = Min., ymax = Max.), fill = colour_line, alpha = alpha) +
  geom_point(data = data_cases_symp, aes(y = value), size = point_size, colour = "black") +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40), expand = c(0, 0)) +
  labs(title = "Symptomatic cases",
       x = "Onset date",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        #strip.background = element_rect(fill = "white"),
        #strip.text = element_text(face = "bold", size = 0),
        axis.title.x.bottom = element_text(vjust = -4, margin = margin(b = 40), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.title.y.right = element_text(margin = margin(l = 20)),
        axis.text.x.bottom = element_text(vjust = -0.4, size = axis_text_size),
        axis.text.y.left = element_text(size = axis_text_size),
        plot.title = element_text(size = title_size, vjust = 5, hjust = title_align, margin = margin(t = 30), face = title_face)) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

ggsave(g, file = here::here("plots", paste0(scenario,"_symp.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

h = ggplot(subset(trajectories_obs, var == "Z_n"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin = Min., ymax = Max.), fill = colour_line, alpha = alpha) +
  geom_point(data = data_cases_non_symp, aes(y = value), size = point_size, colour = "black") +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 5), limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Non-symptomatic cases",
       x = "Test date",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        #strip.background = element_rect(fill = "white"),
        #strip.text = element_text(face = "bold", size = 0),
        axis.title.x.bottom = element_text(vjust = -4, margin = margin(b = 40), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.title.y.right = element_text(margin = margin(l = 20)),
        axis.text.x.bottom = element_text(vjust = -0.4, size = axis_text_size),
        axis.text.y.left = element_text(size = axis_text_size),
        plot.title = element_text(size = title_size, vjust = 5, hjust = title_align, margin = margin(t = 30), face = title_face)) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

ggsave(h, file = here::here("plots", paste0(scenario,"_non_symp.png")), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)


