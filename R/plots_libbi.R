library('ggplot2')
library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")
library("GGally")
library('coda')
library('bayesplot')
library("here")
library("patchwork")
library("data.table")

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

dates <- 
  readr::read_csv(file = here::here("data", "dates_day_no_conv.csv")) %>%
  dplyr::select(dates = date)

scenario = "primary"

posterior = read_libbi(here::here("results", paste0("posterior_",scenario,".rds")))

parameters = summary(posterior,quantiles = c(0.025,0.975))
parameters = parameters[,c("var","2.5%","Median","97.5%")]

prop_trans_total = bi_read(posterior, vars = "P_a")
prop_trans_total = as.data.table(prop_trans_total$P_a)
prop_trans_total = prop_trans_total[time==32]

traces = get_traces(posterior, thin = 100)

#trajectories_resample = predict(posterior, start_time=0, end_time=32, output_every=1)
#trajectories = summary(trajectories_resample, type = "state", quantiles =  c(0.025,0.975))
trajectories = summary(posterior, type = "state", quantiles = c(0.025,0.975))

trajectories_non_obs = subset(trajectories, var %in% c("S","E","I_p","I_sk","I_su","I_a","T","C","R_s","R_n","N_tests"))
trajectories_obs = subset(trajectories, var %in% c("Z_sk","Z_n"))
trajectories_abs_trans = subset(trajectories, var %in% c("dE_a","dE_p","dE_s"))
trajectories_prop_trans = subset(trajectories, var %in% c("dP_a","dP_p","dP_s"))
trajectories_prop_trans_cum = subset(trajectories, var %in% c("P_a","P_p","P_s"))

title_size = 18
axis_title_size = 16
x_axis_text_size = 11
y_axis_text_size = 13
title_align = 0.5
title_face = "bold"
line_thickness = 0.7
colour_line = "#CE2931"
point_size = 1.5
alpha = 0.2
patches_height = 1.6*16
height = 16
aspect_ratio = 1.4

#mcmc_pairs(traces,diag_fun = c("dens"), off_diag_args = list(size = 0.5, alpha = 0.2))

#mcmc_trace(traces)

#mcmc_acf(traces)

# ggplot(trajectories_non_obs, aes(x=time)) +
# geom_line(aes(y=Median), colour = colour_line) +
# geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
# facet_wrap(vars(var), scales="free",nrow=3) +
# ylab("Number")

# ggplot(trajectories_abs_trans, aes(x=time)) +
# geom_line(aes(y=Median), colour = colour_line) +
# geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
# facet_wrap(vars(var), scales="fix",nrow=2) +
# ylab("Proportion") 

# ggplot(trajectories_prop_trans, aes(x=time)) +
# geom_line(aes(y=Median), colour = colour_line) +  
# geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
# facet_wrap(vars(var), scales="fix",nrow=2) +
# ylab("Proportion") 

# ggplot(trajectories_prop_trans_cum, aes(x=time)) +
# geom_line(aes(y=Median), colour = colour_line) +
# geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) + 
# facet_wrap(vars(var), scales="fix",nrow=2) +
# ylab("Proportion") 

parameters

prop_trans_total[,quantile(value,c(0.025,0.5,0.975))]

clinical_fit = ggplot(subset(trajectories_obs, var == "Z_sk"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_cases_symp, aes(y = value), size = point_size, colour = "black") +
  scale_x_continuous(breaks = seq(0, 32, 1), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40), expand = c(0, 0)) +
  labs(title = "Pre or subclinical cases",
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
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

subclinical_fit = ggplot(subset(trajectories_obs, var == "Z_n"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_cases_non_symp, aes(y = value), size = point_size, colour = "black") +
  scale_x_continuous(breaks = seq(0, 32, 1), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Clinical cases",
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
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

tests = ggplot(data_tests, aes(x=time)) + 
  geom_bar(stat="identity", aes(y = value), colour = "white", fill = "#24478f") +
  scale_x_continuous(breaks = seq(0, 32, 1), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 700, 100), limits = c(0, 700), expand = c(0, 0)) +
  labs(title = "Symptom agnostic testing",
       x = "Test date",
       y = "Number of tests") +
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
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none")

fitting = (clinical_fit|subclinical_fit)/(tests|tests) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), plot.tag.position = c(0.03,0.945))

ggsave(fitting, file = here::here("plots","primary", paste0(scenario,"_fitting.png")), width = aspect_ratio*patches_height, height = patches_height, units = "cm", dpi = 300)


relative_inf = ggplot(traces, aes(x=theta_a)) +
  geom_histogram(binwidth = 0.02,fill=colour_line,alpha=0.5) +
  scale_x_continuous(breaks = seq(0, 1.0, 0.1), limits = c(0, 1.05), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 250, 50), limits = c(0, 250), expand = c(0, 0), labels = seq(0,250/5000,50/5000)) +
  labs(title = "",
       x = "Relative infectiousness of subclinical disease",
       y = "Probability") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 15, vjust=2),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0, hjust = 0.5, size = x_axis_text_size, angle = 0),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 0), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

contribution_trans = ggplot(prop_trans_total, aes(x=value)) +
  geom_histogram(binwidth = 0.02,fill=colour_line,alpha=0.5) +
  scale_x_continuous(breaks = seq(0, 1.0, 0.1), limits = c(0, 1.05), expand = c(0, 0), labels = seq(0,100,10)) +
  scale_y_continuous(breaks = seq(0, 25000, 5000), limits = c(0, 25000), expand = c(0, 0), labels = seq(0,25000/500000,5000/500000)) +
  labs(title = "",
       x = "% of transmission from subclinical disease",
       y = "Probability") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 15, vjust=2),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0, hjust = 0.5, size = x_axis_text_size, angle = 0),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 0), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

abs_trans = ggplot(trajectories_abs_trans, aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  scale_x_continuous(breaks = seq(0, 32, 1), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 350, 25), limits = c(0, 350), expand = c(0, 0)) +
  labs(title = "",
       x = "Date",
       y = "Infections per day") +
  facet_wrap(vars(var), scales="fix",nrow=1,labeller = labeller(var=c("dE_a"="Subclinical","dE_p"="Preclinical","dE_s"="Clinical"))) +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 15, vjust=2),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size-1, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 0), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
  theme(legend.position = "none") 

# prop_trans = ggplot(trajectories_prop_trans, aes(x=time)) +
#   geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
#   geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
#   scale_x_continuous(breaks = seq(0, 32, 1), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
#   scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0), labels = seq(0,100,10)) +
#   labs(title = "",
#        x = "Date",
#        y = "% of transmission per day") +
#   facet_wrap(vars(var), scales="fix",nrow=1) + #,labeller = labeller(var=c("dP_a"="Subclinical","dP_p"="Preclinical","dP_s"="Clinical"))) +
#   theme(text = element_text(family = "sans"),
#         panel.background = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.line = element_line(),
#         axis.line.x.bottom = element_line(size = line_thickness),
#         axis.line.y.left = element_line(size = line_thickness),
#         axis.ticks = element_line(size = line_thickness),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(face = "bold", size = 0),
#         axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
#         axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
#         axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size-1, angle = 90),
#         axis.text.y.left = element_text(size = y_axis_text_size),
#         panel.spacing.x=unit(5,"mm"), 
#         plot.title = element_text(size = title_size, vjust = 3, hjust = title_align, margin = margin(t = 0), face = title_face)) +
#   #scale_colour_manual(name = NULL, values = c(colour_no_self_cure,colour_line), labels = c("No self-cure", "Self-cure")) +
#   theme(legend.position = "none") 

trans = (relative_inf|contribution_trans)/abs_trans + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), plot.tag.position = c(0.03,0.97))

ggsave(trans, file = here::here("plots","primary", paste0(scenario,"_trans.png")), width = aspect_ratio*patches_height, height = patches_height, units = "cm", dpi = 300)






