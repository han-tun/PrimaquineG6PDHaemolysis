rm(list=ls())
blues = sapply(seq(0.2,.5,by=0.1), function(x) adjustcolor('blue', alpha.f = x))
load('../Data/HaemolysisData.RData')

Rcpp::sourceCpp('ForwardSim.cpp')

# This looks at the individual functions to get an idea of their behaviour and reasonable parameters
##****** The transit time function for reticulocytes moving from bone marrow to circulation
Hbs = seq(0,20,length.out = 100)
tts = array(dim=length(Hbs))
for(i in 1:length(Hbs)){
  tts[i] = (120-compute_transit_time(C_t_minus_1 = Hb_to_NumberCells(Hb = Hbs[i],
                                                                     BloodVolume = 5,
                                                                     MeanCellHb = 30),
                                     Hb_50_circ = 10, 
                                     k = 10^(-6),
                                     BloodVolume = 5,
                                     MeanCellHb = 30))/24
}
plot(Hbs,tts, type='l', ylab = 'Number of days that retics circulate', 
     xlab='Haemoglobin',main = 'transit time function')


##****** The lifespan function for different mg/kg PMQ doses *******
ds = seq(0,45,length.out = 100)/60
ssls = array(dim = length(ds))
for(i in 1:length(ds)){
  ssls[i] = compute_steady_state_life_span(T_E_star = 115*24,
                                           E_max = .4, # fraction decrease
                                           effectivePMQdose = ds[i],
                                           x_50 = 15/60,
                                           PMQ_slope = 2)
}
plot(ds*60, ssls/24, type='l',xlab = 'Effective PMQ dose in 60 kg adult',
     ylab = 'Lifespan of erythrocytes')

##****** The fold change function: how the production of RBCs changes as a function of Hb *******
Hbs = 1:20
rhos = array(dim = length(Hbs))
for(i in 1:length(Hbs)){
  rhos[i] = Fold_Change_Production(rho_max = 8,
                                   Cstar = Hb_to_NumberCells(Hb = 15,
                                                             BloodVolume = 5,
                                                             MeanCellHb = 30),
                                   C_t_minus_1 = Hb_to_NumberCells(Hb = Hbs[i],
                                                                   BloodVolume = 5,
                                                                   MeanCellHb = 30),
                                   C_50_rho = Hb_to_NumberCells(Hb = 10,
                                                                BloodVolume = 5,
                                                                MeanCellHb = 30))
}
plot(Hbs, rhos, type='l',xlab = 'Hb',ylab = 'Fold change in normoblast production',lwd=3)
abline(h=1, v=15)

##****** Testing out the full simulation function *******
drug_regimen = unlist(lapply(c(7.5, 15, 22.5, 30, 0)/60, function(x) rep(x, 5*24)))
out = forward_sim(drug_regimen = as.double(drug_regimen),
                  rho_max = 5,
                  Hb_steady_state = 15,
                  Hb_50_rho = 10,
                  Hb_50_circ = 10,
                  k = 10^(-6),
                  E_max = .6,
                  x_50 = 15/60,
                  T_m = as.integer(7*24),
                  T_E_star = as.integer(115*24),
                  PMQ_slope = 2,
                  MeanCellHb = 30,
                  BloodVolume = 5)

par(las=1, bty='n', mfrow=c(1,2))
plot((1:length(drug_regimen))/24,out$Hb, lwd=2,#ylim=c(30,45)/3,
     type='l',xlab='days',ylab = 'Haemoglobin (g/dL)')
title("Haemoglobin")
plot((1:length(drug_regimen))/24,out$retic_percent, lwd=2,#ylim=c(0,9),
     type='l',xlab='days',ylab = 'retics (%)')
title("Retic count")


