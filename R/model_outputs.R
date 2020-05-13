library('ggplot2')
library("here")
library("tidyverse")
library("rbi")
library("rbi.helpers")
library('bayesplot')
library("patchwork")
library("data.table")
library("reshape2")
library("parallel")



### Inputs ###

# Posterior to sample from 

posterior_name = "primary"

# LibBi model file to use for sampling 

sample_model_name = "primary_sampler"

# Number of samples 

samples = 100000



### Sample the posterior ###

# Load the model and posterior

model <- rbi::bi_model(file = here::here("bi", paste0(sample_model_name, ".bi")))

posterior = read_libbi(here::here("posteriors", paste0(posterior_name,".rds")))

# Shuffle and thin the posterior to find parameters to sample with

traces_full = data.table(get_traces(posterior))
traces_shuffled = traces_full[c(sample(dim(traces_full)[1],dim(traces_full)[1],replace=FALSE))]
traces_shuffled = traces_shuffled[c(sample(dim(traces_full)[1],dim(traces_full)[1],replace=FALSE))]
traces_shuffled = traces_shuffled[c(sample(dim(traces_full)[1],dim(traces_full)[1],replace=FALSE))]
traces_thinned = traces_shuffled[seq(1,dim(traces_full)[1],dim(traces_full)[1]/samples)]

df_beta_bar =  data.frame(np = seq(1,samples), value = traces_thinned$beta_bar)
df_b_1 =  data.frame(np = seq(1,samples), value = traces_thinned$b_1)
df_tau_11 =  data.frame(np = seq(1,samples), value = traces_thinned$tau_11)
df_tau_22 =  data.frame(np = seq(1,samples), value = traces_thinned$tau_22)
df_chi =  data.frame(np = seq(1,samples), value = traces_thinned$chi)
df_theta_a =  data.frame(np = seq(1,samples), value = traces_thinned$theta_a)
df_theta_p =  data.frame(np = seq(1,samples), value = traces_thinned$theta_p)
df_c_22 =  data.frame(np = seq(1,samples), value = traces_thinned$c_22)

traces_thinned = list(beta_bar = df_beta_bar,
                      b_1 = df_b_1,
                      tau_11 = df_tau_11,
                      tau_22 = df_tau_22,
                      chi = df_chi,
                      theta_a = df_theta_a,
                      theta_p = df_theta_p,
                      c_22 = df_c_22)

remove(traces_shuffled)

# Run the model with the sampled parameters 

trajectories_full = predict(posterior, model = model, init = traces_thinned, nsamples = samples, start_time=0, end_time=32, output_every=1, with=c("transform-obs-to-state"))

remove(traces_thinned)



### Outputs for Figure 1 in main paper ###

trajectories = summary(trajectories_full, type = "state", quantiles = c(0.025,0.975))

prevalence = subset(trajectories, var %in% c("Z_prev"))

observations = summary(trajectories_full, type = c("obs"), quantiles =  c(0.025,0.975))

# Calculate R_0 

NGM_comps_f = c("f_15","f_16","f_17","f_18", "f_19", "f_110", "f_111", "f_112", "f_25", "f_26", "f_27", "f_28", "f_29", "f_210","f_211","f_212")
NGM_comps_v = c("v_11","v_22","v_31","v_33","v_42","v_44","v_53","v_55","v_64","v_66","v_73","v_77","v_84","v_88","v_97","v_99","v_108", "v_1010","v_117","v_1111","v_128","v_1212")
NGM_comp = bi_read(trajectories_full, var = c(NGM_comps_f,NGM_comps_v))
NGM_comp = melt(NGM_comp, id.vars = c("np","time"), variable.name = "L1")
NGM_comp = as.data.table(NGM_comp)
colnames(NGM_comp)=c("np","time","comp","value")
NGM_comp = dcast(NGM_comp, np + time ~ comp)
NGM_comp = as.data.table(NGM_comp)
NGM_comp = NGM_comp[time !=0,]
NGM_comp$R_0 = 0
NGM_comp$number = seq(1,dim(NGM_comp)[1])
R_0_func = function(row){
  f = rep(0,12^2)
  v = f
  f[c(seq(5,12),seq(17,24))] = row[NGM_comps_f]
  f = matrix(f, nrow=12, byrow=TRUE)
  v[c(1,14,25,27,38,40,51,53,64,66,75,79,88,92,103,105,116,118,127,131,140,144)] = row[NGM_comps_v]
  v = matrix(v, nrow=12, byrow=TRUE)
  NGM = f %*% solve(v)
  #print(row["number"])
  return(max(eigen(NGM)$values))
}

cl <- makeCluster(detectCores())
clusterExport(cl, "NGM_comp")
clusterExport(cl, "R_0_func")
clusterExport(cl, "NGM_comps_f")
clusterExport(cl, "NGM_comps_v")
NGM_comp$R_0 = parApply(cl,NGM_comp,1,R_0_func)
stopCluster(cl)

R_0 = NGM_comp[,.("2.5%"= quantile(R_0,0.025), "Median" = quantile(R_0,0.5), "97.5%" = quantile(R_0,0.975)),time]
R_0 = add_row(R_0, .before = 1)
R_0[1,c(2,3,4)] = R_0[2,c(2,3,4)]
R_0[1,1] = 0



### Outputs for Figure 2 in main paper ###

traces = data.table(get_traces(trajectories_full))

proportion_transmission = bi_read(trajectories_full, vars = "P_a")
proportion_transmission = as.data.table(proportion_transmission$P_a)
proportion_transmission = proportion_transmission[time==32]

missed = subset(trajectories, var %in% c("num_s","num_p_and_a") & time == 32)
found = subset(trajectories, var %in% c("R_s","R_n") & time == 32)
missed[,"var"] = c("subclinical","clinical")
found[,"var"] = c("clinical","subclinical")



###  Outputs for text in the main paper ###

parameters = summary(trajectories_full, quantiles = c(0.025,0.975))

chi_numbers = subset(parameters, var == "chi")[,c("2.5%","Median","97.5%")]

prop_transmission_numbers_subclinical = subset(trajectories, var == "P_a" & time == 32)[,c("2.5%","Median","97.5%")]
prop_transmission_numbers_subclinical_IQR = quantile(proportion_transmission[time==32]$value,c(0.25,0.5,0.75))
proportion_transmission_over_half = dim(subset(proportion_transmission, value >= 0.5))[1]/dim(subset(proportion_transmission))[1]

prop_transmission_numbers_preclinical = subset(trajectories, var == "P_p" & time == 32)[,c("2.5%","Median","97.5%")]
prop_transmission_numbers_clinical = subset(trajectories, var == "P_s" & time == 32)[,c("2.5%","Median","97.5%")]

case_infection_ratio = 1/(1-chi_numbers)

infections = subset(trajectories, var %in% c("I_total","per_inf","per_found","per_not_test","per_rec","per_exp","per_found_p_and_a","per_found_s") & time == 32)
infections = infections[,c("var","2.5%","Median","97.5%")]

DIC = DIC(posterior)

R_O_initial = R_0[time==0]
R_O_initial = R_O_initial[,c(2,3,4)]



### Outputs for Table 2 in the main paper

thetas = data.table(thetas = traces$theta_a)

prop_trans_table = proportion_transmission[,.(prop_trans = value)]

R_0_0 = NGM_comp[time == 1, .(R_0_0 = R_0)]

# Calculate R_net for presymptomatic/symptomatic

NGM_comp = bi_read(trajectories_full, var = c(NGM_comps_f,NGM_comps_v))
NGM_comp = melt(NGM_comp, id.vars = c("np","time"), variable.name = "L1")
NGM_comp = as.data.table(NGM_comp)
NGM_comp = NGM_comp[time==1,]
colnames(NGM_comp)=c("np","time","comp","value")
NGM_comp = dcast(NGM_comp, np + time ~ comp)
NGM_comp = as.data.table(NGM_comp)
NGM_comp$R_p = 0
NGM_comp$number = seq(1,dim(NGM_comp)[1])

R_p_func = function(row){
  f = rep(0,12^2)
  v = f
  f[c(seq(5,12),seq(17,24))] = row[NGM_comps_f]
  f = matrix(f, nrow=12, byrow=TRUE)
  v[c(1,14,25,27,38,40,51,53,64,66,75,79,88,92,103,105,116,118,127,131,140,144)] = row[NGM_comps_v]
  v = matrix(v, nrow=12, byrow=TRUE)
  NGM = f %*% solve(v)
  R_p = NGM[1,8]+NGM[2,8]
  print(row["number"])
  return(R_p)
}

cl <- makeCluster(detectCores())
clusterExport(cl, "NGM_comp")
clusterExport(cl, "R_p_func")
clusterExport(cl, "NGM_comps_f")
clusterExport(cl, "NGM_comps_v")
NGM_comp$R_p = parApply(cl,NGM_comp,1,R_p_func)
stopCluster(cl)

R_p = NGM_comp[,.(R_p = R_p)]

theta_table = data.table(thetas,prop_trans_table,R_0_0,R_p) 

theta_results_1 = theta_table[thetas>0 & thetas<=0.01,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                        "prop_trans_u" = quantile(prop_trans,0.975),
                                                        "R_0_0_l" = quantile(R_0_0,0.025),
                                                        "R_0_0_u" = quantile(R_0_0,0.975),
                                                        "R_p_l" = quantile(R_p,0.025),
                                                        "R_p_u" = quantile(R_p,0.975))]

theta_results_2 = theta_table[thetas>0.01 & thetas<=0.25,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                           "prop_trans_u" = quantile(prop_trans,0.975),
                                                           "R_0_0_l" = quantile(R_0_0,0.025),
                                                           "R_0_0_u" = quantile(R_0_0,0.975),
                                                           "R_p_l" = quantile(R_p,0.025),
                                                           "R_p_u" = quantile(R_p,0.975))]

theta_results_3 = theta_table[thetas>0.25 & thetas<=0.5,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                          "prop_trans_u" = quantile(prop_trans,0.975),
                                                          "R_0_0_l" = quantile(R_0_0,0.025),
                                                          "R_0_0_u" = quantile(R_0_0,0.975),
                                                          "R_p_l" = quantile(R_p,0.025),
                                                          "R_p_u" = quantile(R_p,0.975))]

theta_results_4 = theta_table[thetas>0.5 & thetas<=0.75,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                          "prop_trans_u" = quantile(prop_trans,0.975),
                                                          "R_0_0_l" = quantile(R_0_0,0.025),
                                                          "R_0_0_u" = quantile(R_0_0,0.975),
                                                          "R_p_l" = quantile(R_p,0.025),
                                                          "R_p_u" = quantile(R_p,0.975))]

theta_results_5 = theta_table[thetas>0.75 & thetas<=0.99,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                           "prop_trans_u" = quantile(prop_trans,0.975),
                                                           "R_0_0_l" = quantile(R_0_0,0.025),
                                                           "R_0_0_u" = quantile(R_0_0,0.975),
                                                           "R_p_l" = quantile(R_p,0.025),
                                                           "R_p_u" = quantile(R_p,0.975))]

theta_results_6 = theta_table[thetas>0.99 & thetas<=1,.("prop_trans_l" = quantile(prop_trans,0.025),
                                                        "prop_trans_u" = quantile(prop_trans,0.975),
                                                        "R_0_0_l" = quantile(R_0_0,0.025),
                                                        "R_0_0_u" = quantile(R_0_0,0.975),
                                                        "R_p_l" = quantile(R_p,0.025),
                                                        "R_p_u" = quantile(R_p,0.975))]

theta_results_all = rbind(theta_results_1,
                          theta_results_2,
                          theta_results_3,
                          theta_results_4,
                          theta_results_5,
                          theta_results_6)



### Outputs for supplementary materials ###

parameter_results = parameters[,c("var","2.5%","Median","97.5%")]

traces_pairs_1 = traces[seq(1,samples,samples/10000)]
traces_pairs_2 = traces_pairs_1
traces_pairs_1$chain = 1
traces_pairs_2$chain = 2
traces_pairs = rbind(traces_pairs_1,traces_pairs_2)

trans_theta_1 = data.table("Proportion of transmission from asymptomatics" = proportion_transmission$value, "Relative infectiousness of asymptomatics" = traces$theta_a)
trans_theta_2 = trans_theta_1
trans_theta_1$chain = 1
trans_theta_2$chain = 2
trans_theta = rbind(trans_theta_1,trans_theta_2)