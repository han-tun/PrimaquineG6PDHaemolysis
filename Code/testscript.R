rm(list=ls())
blues = sapply(seq(0.2,.5,by=0.1), function(x) adjustcolor('blue', alpha.f = x))
load('../Data/HaemolysisData.RData')

Rcpp::sourceCpp('ForwardSim.cpp')


par(las=1, bty='n')

##** Parameters for the simulations ####
T_E_star = as.integer(115*24)
BloodVolume = 5
MCH = 30
Hb_50_circ = 10
Hb_50_rho = 10
E_max = .3
x_50 = 15/60
PMQ_slope = 2
rho_max = 5
Hb_star = 15


# This looks at the individual functions to get an idea of their behaviour and reasonable parameters
##****** The transit time function for reticulocytes moving from bone marrow to circulation
Hbs = seq(5,16,length.out = 100)
tts = array(dim=length(Hbs))
for(i in 1:length(Hbs)){
  tts[i] = (120-compute_transit_time(C_t_minus_1 = Hb_to_NumberCells(Hb = Hbs[i],
                                                                     BloodVolume = BloodVolume,
                                                                     MeanCellHb = MCH),
                                     Hb_50_circ = Hb_50_circ, 
                                     k = 10^(-6),
                                     BloodVolume = BloodVolume,
                                     MeanCellHb = MCH))/24
}
plot(Hbs,tts, type='l', ylab = 'Number of days that retics circulate', 
     xlab='Haemoglobin',main = 'transit time function')


par(las=1, bty='n', mfrow=c(2,2))

##****** The lifespan function for different mg/kg PMQ doses *******
ds = seq(0,45,length.out = 100)/60
ssls = array(dim = length(ds))
for(i in 1:length(ds)){
  ssls[i] = compute_steady_state_life_span(T_E_star = T_E_star,
                                           E_max = E_max, # fraction decrease
                                           effectivePMQdose = ds[i],
                                           x_50 = x_50,
                                           PMQ_slope = PMQ_slope)
}
plot(ds*60, ssls/24, type='l',xlab = 'Effective PMQ dose in 60 kg adult',
     ylab = 'Lifespan of erythrocytes')

##****** The fold change function: how the production of RBCs changes as a function of Hb *******
Hbs = 1:20
rhos = array(dim = length(Hbs))
for(i in 1:length(Hbs)){
  rhos[i] = Fold_Change_Production(rho_max = rho_max,
                                   Cstar = Hb_to_NumberCells(Hb = Hb_star,
                                                             BloodVolume = BloodVolume,
                                                             MeanCellHb = MCH),
                                   C_t_minus_1 = Hb_to_NumberCells(Hb = Hbs[i],
                                                                   BloodVolume = BloodVolume,
                                                                   MeanCellHb = MCH),
                                   C_50_rho = Hb_to_NumberCells(Hb = Hb_50_rho,
                                                                BloodVolume = BloodVolume,
                                                                MeanCellHb = MCH))
}
plot(Hbs, rhos, type='l',xlab = 'Hb',ylab = 'Fold change in normoblast production',lwd=3)
abline(h=1, v=15)

##****** Testing out the full simulation function *******
drug_regimen = unlist(lapply(c(7.5, 15, 22.5, 30, 0, 0)/60, function(x) rep(x, 5*24)))
out = forward_sim(drug_regimen = as.double(drug_regimen),
                  rho_max = rho_max,
                  Hb_steady_state = Hb_star,
                  Hb_50_rho = Hb_50_rho,
                  Hb_50_circ = Hb_50_circ,
                  k = 10^(-6),
                  E_max = E_max,
                  x_50 = x_50,
                  T_m = as.integer(7*24),
                  T_E_star = T_E_star,
                  PMQ_slope = PMQ_slope,
                  MeanCellHb = MCH,
                  BloodVolume = BloodVolume)

plot((1:length(drug_regimen))/24,out$Hb, lwd=2,
     ylim=range(c(out$Hb, Haemodata_Analysis$Hb)),
     type='l',xlab='days',ylab = 'Haemoglobin (g/dL)')
title("Haemoglobin")

mean_Hbs = array(dim=length(unique(Haemodata_Analysis$Day)))
mean_Retics = array(dim=length(unique(Haemodata_Analysis$Day)))
for(dd in unique(Haemodata_Analysis$Day)){
  ind = Haemodata_Analysis$Day==dd
  points(dd, mean(Haemodata_Analysis$Hb[ind]), pch=18, col='red')
  lines(c(dd,dd), quantile(Haemodata_Analysis$Hb[ind], probs = c(.1,.9)))
  mean_Hbs[which(dd==unique(Haemodata_Analysis$Day))] = mean(Haemodata_Analysis$Hb[ind])
}
lines(unique(Haemodata_Analysis$Day), mean_Hbs, col='red', lwd=2)


plot((1:length(drug_regimen))/24,out$retic_percent, lwd=2,
     ylim=range(c(Haemodata_Analysis$Rectic_percent, out$retic_percent), na.rm = T),
     type='l',xlab='days',ylab = 'retics (%)')
for(dd in unique(Haemodata_Analysis$Day)){
  ind = Haemodata_Analysis$Day==dd
  points(dd, mean(Haemodata_Analysis$Rectic_percent[ind]), pch=18, col='red')
  lines(c(dd,dd), quantile(Haemodata_Analysis$Rectic_percent[ind], probs = c(.1,.9),na.rm = T))
  mean_Retics[which(dd==unique(Haemodata_Analysis$Day))] = mean(Haemodata_Analysis$Rectic_percent[ind])
}
title("Retic count")


