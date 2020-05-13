library("here")

### Data preparation ###

# Symptomatic data 

data_symp = read.csv(file = here::here("data", "raw_data_symp_nishiura_20.csv"), header = TRUE, sep = ",")    
data_symp$all = rowSums(data_symp[,c(-1)])                                                
data_symp$pass = data_symp$pass_cc + data_symp$pass_other
data_symp$day_no = seq(1,length(data_symp$onset_date),1)                                  

data_cases_symp = data_symp[,c("onset_date","day_no","all","pass","crew")]     
write.csv(data_cases_symp, file = here::here("data", "data_cases_symp.csv"), row.names = FALSE)         

# Non-symptomatic data

data_non_symp = read.csv(file = here::here("data", "raw_data_non_symp_mizumoto_20.csv"), header = TRUE, sep = ",")
data_non_symp$tests_non_symp = data_non_symp$tests - data_non_symp$symp                       
data_non_symp$day_no = seq(1,length(data_symp$onset_date),1)                                  
data_non_symp$tests_non_symp_cum = cumsum(data_non_symp$tests_non_symp)                       

data_cases_non_symp = data_non_symp[,c("test_date","day_no","non_symp")]     
write.csv(data_cases_non_symp, file = here::here("data", "data_cases_non_symp.csv"), row.names = FALSE)      

# Test data

data_tests = data_non_symp[,c("test_date","day_no","tests_non_symp","tests_non_symp_cum")]    
write.csv(data_tests, file = here::here("data", "data_tests.csv"), row.names = FALSE)

       

