setwd("~/Dropbox/MORU/RBC modelling/PrimaquineHaemolysis/PrimaquineG6PDHaemolysis/Code")
library(dplyr)
library(tictoc)

rm(list=ls())
load('../../HaemolysisData.RData')
Rcpp::sourceCpp('ForwardSim.cpp')
source('inference_functions.R')
source('optimisation_functions.R')
par(las=1, bty='n')


##** Parameters for the simulations ####
params = unlist(list(T_E_star = as.integer(115*24),
                     BloodVolume = 5,
                     MCH = 30,
                     Hb_50_circ = 10,
                     Hb_50_rho = 10,
                     E_max = .3,
                     x_50 = 15/60,
                     PMQ_slope = 3,
                     rho_max = 5,
                     Hb_star = 14.5))

Haemodata_Analysis = filter(Haemodata_Analysis, Day >= 0 )
Haemodata_Analysis$Day = Haemodata_Analysis$Day + 1


#*********** Run accept-reject style simulation by generating from the prior *************
drug_regimen = array(sapply(c(7.5,15,22.5,30,0,0),rep,24*5))/60

N = 100000
x = array(dim = N)
params_list = generate_from_prior(n = N)
tic()
pb=txtProgressBar(max = N, style = 3)
for(i in 1:N){
  x[i] = compute_post_density(params = params_list[[i]],drug_regimen = drug_regimen,
                              Haemodata_Analysis = Haemodata_Analysis)
  setTxtProgressBar(pb, value = i)
}
toc()

hist(x[x > -1000])
params_max = params_list[[which.max(x)]]
out_max = estimate_sim_profile(params = params_max,drug_regimen = drug_regimen)


par(mfrow=c(1,2),las=1,bty='n')
sim_black = adjustcolor('black',alpha.f = .1)
plot((1:length(drug_regimen))/24,out_max$Hb, lwd=4,
     ylim=range(c(out_max$Hb, Haemodata_Analysis$Hb)),
     type='l',xlab='days',ylab = 'Haemoglobin (g/dL)')
plot_Hb_data(Haemodata_Analysis, col_line = 'red')

title("Haemoglobin")
ind = sample(which(x > quantile(x, probs = 0.999)),replace = T,size = 200)
for(j in ind){
  params = params_list[[j]]
  out = estimate_sim_profile(params = params,drug_regimen = drug_regimen)
  lines((1:length(drug_regimen))/24,out$Hb, col = sim_black)
}


plot((1:length(drug_regimen))/24,out$retic_percent, lwd=4,
     ylim=range(c(Haemodata_Analysis$Rectic_percent, out$retic_percent), na.rm = T),
     type='l',xlab='days',ylab = 'retics (%)')
ind = sample(which(x > quantile(x, probs = 0.999)),replace = T,size = 200)
plot_Retic_data(Haemodata_Analysis,col_line = 'red')
for(j in ind){
  params = params_list[[j]]
  out = estimate_sim_profile(params = params,drug_regimen = drug_regimen)
  lines((1:length(drug_regimen))/24,out$retic_percent, col = sim_black)
}
title("Reticulocyte count")

par(mfrow=c(1,1))
plot_dose_response(params_max,lwd=4, ylim = c(50,120),add = F)
for(j in ind){
  params = params_list[[j]]
  plot_dose_response(params,lwd=1, col=sim_black,add=T)
}

optimal_two_week_reg = optimise_regimen(params = params_max,min_total_dose = 375,duration = 14)
optimal_20day_reg = optimise_regimen(params = params_max,min_total_dose = 375,duration = 20)

par(mfrow=c(1,2))

drug_regimen = array(sapply(c(7.5,15,22.5,30,0,0),rep,24*5))/60
drug_regimen_14days = array(sapply(c(optimal_two_week_reg,rep(0,12)), rep, 24))/60
drug_regimen_20days = array(sapply(c(optimal_20day_reg,rep(0,5)), rep, 24))/60

# Plot the regimens
ind=drug_regimen_20days>0
plot(1:(24*20), drug_regimen_20days[ind]*60, type = 'l', lwd=3, xaxt='n',ylab='Daily mg dose',xlab='Days',
     ylim=c(5,45))
axis(1, at = c(0,5,10,15,20)*24,labels = c(1,5,10,15,20))
ind=drug_regimen_14days>0
lines(1:(24*14), drug_regimen_14days[ind]*60, lwd=3, col='blue')
ind=drug_regimen>0
lines(1:(24*20), drug_regimen[ind]*60, lwd=3, lty=2)
legend('bottomright',col = c('black','black','blue'), 
       legend = c('Currently tested', 'Predicted optimal: 20days','Predicted optimal: 14days'),
       lwd=3, lty = c(2,1,1), inset=0.01)

out_14days=estimate_sim_profile(params = params_max,drug_regimen = drug_regimen_14days)
out_20days=estimate_sim_profile(params = params_max,drug_regimen = drug_regimen_20days)
out =estimate_sim_profile(params = params_max,drug_regimen = drug_regimen)

plot((1:length(drug_regimen_14days))/24,out_14days$Hb, lwd=4,
     ylim=range(c(out_14days$Hb, Haemodata_Analysis$Hb)),
     type='l',xlab='Days',ylab = 'Haemoglobin (g/dL)',col='blue')
lines((1:length(drug_regimen_20days))/24,out_20days$Hb, lwd=4)
lines((1:length(drug_regimen))/24,out$Hb, lwd=4, lty=2)


##********** Primaquine followed by tafenoquine ***************####
optimal_7day_reg = optimise_regimen(params = params_max,min_total_dose = 50,duration = 7, min_increment = 2.5, max_single_dose=15)
optimal_5day_reg = optimise_regimen(params = params_max,min_total_dose = 50,duration = 5, min_increment = 2.5, max_single_dose=15)
wait_one_week = rep(0,7)
tafenoquine_350 = c(rep(20,3*7), rep(0,10)) # equivalence assumed....


ll=1
par(mfrow=c(1,2),las=1)
plot(1:(24*7), array(sapply(optimal_7day_reg,rep,24)), yaxt='n', type='l',xaxt='n',
     ylab='Primaquine dose', xlab='Days', lwd=3, col=adjustcolor('red',ll))
axis(1, at = (0:7)*24, labels = 0:7)
axis(2, at =c(5,7.5,10,15))
lines(1:(24*5), array(sapply(optimal_5day_reg,rep,24)), col = adjustcolor('green',ll),lwd=3)

comb = list(TQ = tafenoquine_350,
            comb_1 =c(optimal_7day_reg, wait_one_week, tafenoquine_350),
            comb_2 = c(optimal_5day_reg, wait_one_week, tafenoquine_350),
            comb_3 = c(optimal_7day_reg, tafenoquine_350),
            comb_4 = c(optimal_5day_reg, tafenoquine_350))
cols = c('blue', rep(c(adjustcolor('red',ll), adjustcolor('green',ll)),2))
ltys = c(1, 1,1,2,2)
plot(NA,NA, xlim = c(0,30), ylim = c(12,15), xlab = 'Days from start of regimen',
     ylab = 'Haemoglobin (g/dL)')
for(i in 1:length(comb)){
  drug_regimen = array(sapply(comb[[i]],rep,24))/60
  out = estimate_sim_profile(params = params_max,drug_regimen = drug_regimen)
  lines((1:length(drug_regimen))/24,out$Hb,lwd=3, col=cols[i],lty=ltys[i])
}
legend('topright', col = cols, lty=ltys, lwd=3,
       legend = c('Tafenoquine alone',
                  '7 days PMQ + TQ (week wait)',
                  '5 days PMQ + TQ (week wait)',
                  '7 days PMQ + TQ (no wait)',
                  '5 days PMQ + TQ (no wait)'))

par(las=1, mfrow=c(1,1))
params_illsu = params_max
params_illsu['T_E_star'] = 115*24
plot_dose_response(params_illsu,lwd=3)


