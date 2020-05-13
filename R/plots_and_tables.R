library('ggplot2')
library("here")
library("tidyverse")
library('bayesplot')
library("patchwork")
library("data.table")
library("reshape2")

# Run 'model_outputs.R' in the same session and set up folder 'results' prior to running this script 



### Results in the main paper ###

# Results in the text

num_res = 17

text_results = data.frame("result" = rep(0,num_res), "lower" = rep(0,num_res),"median" = rep(0,num_res), "upper" = rep(0,num_res), "num" = rep(0,num_res))

text_results[1,] = c("chi",100*round(chi_numbers,2),0)
text_results[2,] = c("prop_trans_sub",100*round(prop_transmission_numbers_subclinical,2),0)
text_results[3,] = c("prop_trans_sub_over_half",rep(0,3),100*round(proportion_transmission_over_half,2))
text_results[4,] = c("prop_trans_pre",100*round(prop_transmission_numbers_preclinical,2),0)
text_results[5,] = c("prop_trans_clin",100*round(prop_transmission_numbers_clinical,2),0)
text_results[6,] = c("case_inf_ratio",round(case_infection_ratio,1),0)
text_results[seq(7,14),1] = infections$var
text_results[7,seq(2,5)] = data.frame(round(infections[1,seq(2,4)],0),rep(0,1))
text_results[seq(8,14),seq(2,5)] = data.frame(100*round(infections[seq(2,8),seq(2,4)],2),rep(0,7))
text_results[15,] = c("DIC",rep(0,3),round(DIC,0))
text_results[16,] = c("prop_trans_sub_IQR",100*round(prop_transmission_numbers_subclinical_IQR,2),0)
text_results[17,] = c("R_O",round(R_O_initial,1),0)

write.csv(text_results, file = here::here("results", paste0(posterior_name), "text_results.csv"), row.names = FALSE) 



# Table 2 

theta_results = data.table("theta_int" = rep("blank",6),theta_results_all)

theta_results[1,1] = "0-1%"
theta_results[2,1] = "1-25%"
theta_results[3,1] = "25-50%"
theta_results[4,1] = "50-75%"
theta_results[5,1] = "75-99%"
theta_results[6,1] = "99-100%"

theta_results[,c(2,3)] = 100*round(theta_results[,c(2,3)],2)
theta_results[,c(4,5,6,7)] = round(theta_results[,c(4,5,6,7)],1)

write.csv(theta_results, file = here::here("results", paste0(posterior_name), "table.csv"), row.names = FALSE) 



### Data for figures ###

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

# Dates 

dates <- 
  readr::read_csv(file = here::here("data", "dates_day_no_conv.csv")) %>%
  dplyr::select(dates = date)

dates = dates[seq(1,33,2),]

# Key dates 

index_case_pos = 13
quarantine_starts = 16

# Non-symptomatic data in prevalence format

data_prev = data_cases_non_symp
data_prev$value = data_prev$value/data_tests$value
data_prev$lower_unc = qbinom(c(0.025),data_tests$value,data_prev$value)/data_tests$value
data_prev$upper_unc = qbinom(c(0.975),data_tests$value,data_prev$value)/data_tests$value



### Figure settings ###

title_size = 13
axis_title_size = 16
x_axis_text_size = 12
y_axis_text_size = 13
title_align = 0.5
title_face = "bold"
line_thickness = 0.6
colour_line = "#CE2931"
point_size = 1
alpha = 0.2
patches_height = 1.6*16
height = 16
aspect_ratio = 1.4
error_bar_width = 1.4
date_line_thickness = 0.4
date_line_colour = "gray35"
prior_colour = "deepskyblue2"



### Figures in the main paper ###

# Figure 1

clinical_fit_1 = 
  ggplot(subset(observations, var == "symp_1"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_cases_symp_1, aes(y = value), size = point_size, colour = "black") +
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40), expand = c(0, 0)) +
  labs(title = "Confirmed symptomatic cases (Crew)",
       x = "Date of symptom onset",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

clinical_fit_2 = 
  ggplot(subset(observations, var == "symp_2"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_cases_symp_2, aes(y = value), size = point_size, colour = "black") +
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(0, 40), expand = c(0, 0)) +
  labs(title = "Confirmed symptomatic cases (Passengers)",
       x = "Date of symptom onset",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

subclinical_fit = 
  ggplot(subset(observations, var == "non_symp"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_cases_non_symp, aes(y = value), size = point_size, colour = "black") +
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  labs(title = "Confirmed pre/asymptomatic cases (All)",
       x = "Date of test",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none") 

prevalence_plot =
  ggplot(subset(prevalence, var == "Z_prev"), aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_point(data = data_prev, aes(y = value), size = point_size, colour = "black") +
  geom_errorbar(data = data_prev, aes(ymin = lower_unc, ymax = upper_unc), size = line_thickness-0.2, width = 0.5) + 
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.05), limits = c(0, 0.4), expand = c(0, 0)) +
  labs(title = "Prevalence of pre/asymptomatic individuals",
       x = "Date of test",
       y = "Proportion") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none") 

basic_reproduction_number =
  ggplot(R_0, aes(x=time)) +
  geom_line(aes(y = Median), size = line_thickness, na.rm=TRUE, colour = colour_line) +
  geom_ribbon(aes(ymin=`2.5%`, ymax= `97.5%`), fill = colour_line, alpha=alpha) +
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 30, 2), limits = c(0, 30), expand = c(0, 0)) +
  labs(title = "Basic reproduction number: All",
       x = "Date",
       y = "Basic reproduction number") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

tests = ggplot(data_tests, aes(x=time)) + 
  geom_bar(stat="identity", aes(y = value), colour = "white", fill = "#24478f") +
  scale_x_continuous(breaks = seq(0, 32, 2), limits = c(0, 32.5), expand = c(0, 0), labels = dates$dates) +
  scale_y_continuous(breaks = seq(0, 700, 100), limits = c(0, 700), expand = c(0, 0)) +
  geom_vline(xintercept = index_case_pos, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  geom_vline(xintercept = quarantine_starts, linetype="dashed", color = date_line_colour, size = date_line_thickness) +
  labs(title = "Symptom agnostic testing",
       x = "Date of test",
       y = "Number of tests") +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(),
        axis.line.x.bottom = element_line(size = line_thickness),
        axis.line.y.left = element_line(size = line_thickness),
        axis.ticks = element_line(size = line_thickness),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 20), size = axis_title_size),
        axis.title.y.left = element_text(vjust= 3, margin = margin(l = 15), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.3, hjust = 0,size = x_axis_text_size, angle = 90),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size, vjust = 6, hjust = title_align, margin = margin(t = 20), face = title_face)) +
  theme(legend.position = "none")

figure_1 = (clinical_fit_1|clinical_fit_2|subclinical_fit) / (prevalence_plot|basic_reproduction_number|tests) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), plot.tag.position = c(0.03,0.97))

ggsave(figure_1, file = here::here("results", paste0(posterior_name), "figure_1.png"), width = aspect_ratio*patches_height, height = patches_height, units = "cm", dpi = 300)

# Figure 2 

prior_chi = data.frame(chi=seq(0,samples)/(1+samples))

proportion_subclinical_distribution = 
  ggplot() +
  geom_density(data = prior_chi, aes(x=chi), fill=prior_colour, alpha=0.4, colour = NA) +
  geom_density(data = traces, aes(x=chi), fill=colour_line, alpha=0.6, colour = NA) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0),labels = seq(0,1,0.1)) +
  scale_y_continuous(breaks = seq(0, 24, 6), limits = c(0, 24), expand = c(0, 0)) +
  labs(title = "Proportion of infections that are asymptomatic",
       x = "Proportion",
       y = "Density") +
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
        plot.title = element_text(size = title_size+2, vjust = 6, hjust = title_align, margin = margin(t = 10), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  theme(legend.position = "none") 

prior_theta_a = data.frame(theta_a=seq(0,samples)/(1+samples))

relative_infectiousness_distribution = 
  ggplot() +
  geom_density(data = prior_theta_a, aes(x=theta_a), fill=prior_colour, alpha=0.4, colour = NA) +
  geom_density(data = traces, aes(x=theta_a), fill=colour_line, alpha=0.6, colour = NA) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0),labels = seq(0,1,0.1)) +
  scale_y_continuous(breaks = seq(0, 2.0, 0.5), limits = c(0, 2.0), expand = c(0, 0)) +
  labs(title = "Relative infectiousness of asymptomatics",
       x = "Proportion",
       y = "Density") +
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
        plot.title = element_text(size = title_size+2, vjust = 6, hjust = title_align, margin = margin(t = 10), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  theme(legend.position = "none") 

contribution_to_transmission_distribution = 
  ggplot() +
  geom_density(data = proportion_transmission, aes(x=value), fill=colour_line, alpha=0.6, colour = NA) +
  geom_vline(xintercept = prop_transmission_numbers_subclinical_IQR[[1]], linetype="dashed", color = date_line_colour, size = date_line_thickness+0.4) +
  geom_vline(xintercept = prop_transmission_numbers_subclinical_IQR[[2]], linetype="longdash", color = date_line_colour, size = date_line_thickness+0.4) +
  geom_vline(xintercept = prop_transmission_numbers_subclinical_IQR[[3]], linetype="dashed", color = date_line_colour, size = date_line_thickness+0.4) +
  scale_x_continuous(breaks = seq(0, 1.0, 0.1), limits = c(0, 1.05), expand = c(0, 0), labels = seq(0,1,0.1)) +
  scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
  labs(title = "Proportion of transmission from asymptomatics",
       x = "Proportion",
       y = "Density") +
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
        plot.title = element_text(size = title_size+2, vjust = 6, hjust = title_align, margin = margin(t = 10), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  theme(legend.position = "none") 

cases_missed_and_found =
  ggplot() +
  geom_bar(data = missed, stat="identity", aes(x = var, y = Median), alpha = 0.3, fill = colour_line,  width = 0.75) +
  geom_bar(data = found, stat="identity", aes(x = var, y = Median), alpha = 0.6, fill = colour_line,  width = 0.75) +
  geom_errorbar(data = rbind(missed,found), aes(x = var, ymin = `2.5%`, ymax = `97.5%`), width = 0.1, position = "dodge2") +
  scale_x_discrete(labels = c("Symptomatic","Pre/asymptomatic"), expand = c(0.1, 0.5)) + # breaks = seq(start_year, end_year, 25), limits = c(start_year, end_year)) +
  scale_y_continuous(breaks = seq(0, 1200, 200), limits = c(0, 1200), expand = c(0, 0)) +
  labs(title = "Total and confirmed infections",
       x = "",
       y = "Number") +
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
        axis.text.x.bottom = element_text(vjust = -1.5, hjust = 0.5, size = x_axis_text_size + 5, angle = 0),
        axis.text.y.left = element_text(size = y_axis_text_size),
        plot.title = element_text(size = title_size+2, vjust = 6, hjust = title_align, margin = margin(t = 10), face = title_face),
        panel.spacing.x=unit(5,"mm")) +
  theme(legend.position = "none")

figure_2 = (proportion_subclinical_distribution|relative_infectiousness_distribution)/(cases_missed_and_found|contribution_to_transmission_distribution) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16), plot.tag.position = c(0.03,1))

ggsave(figure_2, file = here::here("results", paste0(posterior_name), "figure_2.png"), width = 0.9*aspect_ratio*patches_height, height = patches_height, units = "cm", dpi = 300)



### Supplementary materials ###

# Parameter estimates - Supplementary Table 3

parameter_results_table = as.data.table(parameter_results)

parameter_results_table[var=="beta_bar",c(2,3,4)] = round(parameter_results_table[var=="beta_bar",c(2,3,4)],2)
parameter_results_table[var=="c_22",c(2,3,4)] = round(parameter_results_table[var=="c_22",c(2,3,4)],2)
parameter_results_table[var=="b_1",c(2,3,4)] = round(parameter_results_table[var=="b_1",c(2,3,4)],2)
parameter_results_table[var=="tau_11",c(2,3,4)] = round(parameter_results_table[var=="tau_11",c(2,3,4)],1)
parameter_results_table[var=="tau_22",c(2,3,4)] = round(parameter_results_table[var=="tau_22",c(2,3,4)],1)
parameter_results_table[var=="chi",c(2,3,4)] = round(parameter_results_table[var=="chi",c(2,3,4)],2)
parameter_results_table[var=="theta_a",c(2,3,4)] = round(parameter_results_table[var=="theta_a",c(2,3,4)],2)
parameter_results_table[var=="theta_p",c(2,3,4)] = round(parameter_results_table[var=="theta_p",c(2,3,4)],2)

write.csv(parameter_results_table, file = here::here("results", paste0(posterior_name), "parameters.csv"), row.names = FALSE) 

# Trace plot - Supplementary Figure 3

bayesplot_theme_set(theme_default())

bayesplot_theme_update(text = element_text(family = "sans"),
                       axis.text.x.bottom = element_text(vjust = 0, hjust = 0.5, size = 9, angle = 0))

trace_plot = 
  mcmc_trace(traces_full) +
  scale_x_continuous(breaks = seq(0, dim(traces_full)[1], dim(traces_full)[1]/4), limits = c(0, dim(traces_full)[1]), labels = scales::scientific) 

ggsave(trace_plot, file = here::here("results", paste0(posterior_name), "trace_plot.png"), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

# Correlation plot - Supplementary Figure 2

bayesplot_theme_set(theme_default())

bayesplot_theme_update(text = element_text(family = "sans"),
                       axis.text.x.bottom = element_text(vjust = 0, hjust = 0.5, size = 9, angle = 0),
                       axis.text.y.left = element_text(vjust = 0, hjust = 0, size = 9, angle = 0))

pairs_plot = mcmc_pairs(traces_pairs,diag_fun = c("hist"), off_diag_args = list(size = 0.5, alpha = 0.035))

ggsave(pairs_plot, file = here::here("results", paste0(posterior_name), "pairs_plot.png"), width = aspect_ratio*1.15*patches_height, height = 1.15*patches_height, units = "cm", dpi = 300)

# Transmission/theta correlation plot - Supplementary Figure 4 

bayesplot_theme_set(theme_default())

bayesplot_theme_update(text = element_text(family = "sans"),
                       axis.text.x.bottom = element_text(vjust = 0, hjust = 0.5, size = 11, angle = 0),
                       axis.text.y.left = element_text(vjust = 0, hjust = 0, size = 11, angle = 0))

corr_plot = 
  mcmc_pairs(trans_theta,diag_fun = c("hist"), off_diag_args = list(size = 0.35, alpha = 0.035)) 

ggsave(corr_plot, file = here::here("results", paste0(posterior_name), "corr_plot.png"), width = aspect_ratio*height, height = height, units = "cm", dpi = 300)



