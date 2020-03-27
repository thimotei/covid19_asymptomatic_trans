setwd("/Users/jonemery/Google Drive/LSHTM - Research Fellow/03 R/Working directory")
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

### Data preparation ###

### Symptomatics 

# Number of symptomatic cases by onset date for passengers (cabin contacts), passengers (other) and crew 
# extracted from figure 1 in https://www.mdpi.com/2077-0383/9/3/657 

data_symp = read.csv(file = "raw_data_symp_nishiura_20.csv", header = TRUE, sep = ",")    # raw data 
data_symp$all = rowSums(data_symp[,c(-1)])                                                # aggregation into total cases
data_symp$day_no = seq(1,length(data_symp$onset_date),1)                                  # add in day numbers (day 1 = Jan 20th)

data_cases_symp = data_symp[,c("onset_date","day_no","all")]     # Incident symptomatic cases by onset date for model fitting
write.csv(data_cases_symp, file = "data_cases_symp.csv", row.names = FALSE)             

### Non-symptomatics 

# Total number of tests, positive tests, symptomatic cases and non-symptomatic cases by date of PCR test 
# extracted from Table 1 in https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180

data_non_symp = read.csv(file = "raw_data_non_symp_mizumoto_20.csv", header = TRUE, sep = ",") # raw data
data_non_symp$tests_non_symp = data_non_symp$tests - data_non_symp$symp                        # number of tests performed on everyone except symptomatics 
data_non_symp$day_no = seq(1,length(data_symp$onset_date),1)                                   # add in day numbers (day 1 = Jan 20th)
data_non_symp$tests_non_symp_cum = cumsum(data_non_symp$tests_non_symp)                        # cumulatiec number of tests performed on everyone except symptomatics 

data_tests = data_non_symp[,c("test_date","day_no","tests_non_symp_cum")]    # Cumulative PCR tests as input for PCR screening function in the model
write.csv(data_tests, file = "data_tests.csv", row.names = FALSE)

data_cases_non_symp = data_non_symp[,c("test_date","day_no","non_symp")]     # Incident non-symptomatic cases by date of PCR test for model fitting
write.csv(data_cases_non_symp, file = "data_cases_non_symp.csv", row.names = FALSE)             


### Plots ###

height = 16
aspect_ratio = 1.4
title_size = 18
title_face = "bold"
axis_title_size = 16
x_axis_text_size = 10
y_axis_text_size = 13
axis_line_thickness = 1.1

### Symptomatics ###

ggplot() +
  geom_bar(data = data_symp, stat="identity", aes(x = day_no, y = all), colour = "white", fill = "#24478f") +
  scale_x_continuous(labels = data_symp$onset_date, breaks = data_symp$day_no, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,34,2), limits = c(0,34), expand = c(0,0)) +
  labs(title = "Symptomatic cases with reported onset date (n=199) ",
       x = "Date of symptom onset",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        plot.title = element_text(size = title_size, vjust = 8, hjust = 0.5, margin = margin(t = 30), face = title_face),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 18), size = axis_title_size),
        axis.title.y.left = element_text(vjust = 4, margin = margin(l = 28), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.5, size = x_axis_text_size, angle = 90, margin = margin(t = 6)),
        axis.text.y.left = element_text(size = y_axis_text_size),
        axis.line.x.bottom = element_line(size = axis_line_thickness),
        axis.line.y.left = element_line(size = axis_line_thickness),
        axis.ticks = element_line(size = axis_line_thickness),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("data_symp_nishiura_20.png", width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

### Non-symptomatics ###

ggplot() +
  geom_bar(data = data_non_symp, stat="identity", aes(x = day_no, y = non_symp), colour = "white", fill = "#24478f") +
  scale_x_continuous(labels = data_non_symp$test_date, breaks = data_non_symp$day_no, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,70,5), limits = c(0,70), expand = c(0,0)) +
  labs(title = "Non-symptomatic cases by test date (n=320)",
       x = "Date of test",
       y = "Number of cases") +
  theme(text = element_text(family = "sans"),
        plot.title = element_text(size = title_size, vjust = 8, hjust = 0.5, margin = margin(t = 30), face = title_face),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 18), size = axis_title_size),
        axis.title.y.left = element_text(vjust = 4, margin = margin(l = 28), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.5, size = x_axis_text_size, angle = 90, margin = margin(t = 6)),
        axis.text.y.left = element_text(size = y_axis_text_size),
        axis.line.x.bottom = element_line(size = axis_line_thickness),
        axis.line.y.left = element_line(size = axis_line_thickness),
        axis.ticks = element_line(size = axis_line_thickness),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("data_non_symp_mizumoto_20.png", width = aspect_ratio*height, height = height, units = "cm", dpi = 300)

### Testing (excluding on symptomatics) ###

ggplot() +
  geom_bar(data = data_non_symp, stat="identity", aes(x = day_no, y = tests_non_symp_cum), colour = "white", fill = "#24478f") +
  scale_x_continuous(labels = data_non_symp$test_date, breaks = data_non_symp$day_no, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,3000,250), limits = c(0,3000), expand = c(0,0)) +
  labs(title = "Cumulative tests (excluding symptomatics) by test date",
       x = "Date of test",
       y = "Cumulative number of tests") +
  theme(text = element_text(family = "sans"),
        plot.title = element_text(size = title_size, vjust = 8, hjust = 0.5, margin = margin(t = 30), face = title_face),
        axis.title.x.bottom = element_text(vjust = -2, margin = margin(b = 18), size = axis_title_size),
        axis.title.y.left = element_text(vjust = 4, margin = margin(l = 28), size = axis_title_size),
        axis.text.x.bottom = element_text(vjust = 0.5, size = x_axis_text_size, angle = 90, margin = margin(t = 6)),
        axis.text.y.left = element_text(size = y_axis_text_size),
        axis.line.x.bottom = element_line(size = axis_line_thickness),
        axis.line.y.left = element_line(size = axis_line_thickness),
        axis.ticks = element_line(size = axis_line_thickness),
        panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("data_tests_mizumoto_20.png", width = aspect_ratio*height, height = height, units = "cm", dpi = 300)





