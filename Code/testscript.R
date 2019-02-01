rm(list=ls())
blues = sapply(seq(0.2,.5,by=0.1), function(x) adjustcolor('blue', alpha.f = x))
load('HaemolysisData.RData')
drug_regimen = unlist(lapply(c(7.5, 15, 22.5, 30, 0), 
                             function(x) rep(x, 5*24)))
days = 1:(length(drug_regimen)/24)
Hbs = Retics = array(dim = length(days))
for(d in days){
  Hbs[d] = mean(Haemodata_Analysis$Hb[Haemodata_Analysis$Day==d])
  Retics[d] = mean(Haemodata_Analysis$Rectic_percent[Haemodata_Analysis$Day==d])
}
Rcpp::sourceCpp('ForwardSim.cpp')

#drug_regimen = rep(25,24*2 + 1)
out=forward_sim(drug_regimen = as.double(drug_regimen),
                max_production = 5,
                mid_haematocrit = 30,
                HCT_steady_state = 45,
                mid_transit = 35,
                hill_coef = 5,
                time_lag = 5*24,
                T_lifeMin = 50*24,
                PMQ_slope = 2, 
                Effect50_dose = 30,
                T_lifeSpan_norm = 100*24)
plot(out$effect, type='l')
plot((1:length(out$CiculatingRBCs))/24, out$CiculatingRBCs/10^5,pch='.')
abline(v= 120*24 + out$T_retic)
par(las=1, bty='n', mfrow=c(1,2))
plot((1:length(drug_regimen))/24,out$haemtocrit, lwd=3,
     type='l',xlab='days',ylab = 'HCT(%)', ylim = c(30,45))
polygon(c(0,5,5,0), c(0,0,1000,1000), col = blues[1], border = NA)
polygon(c(5,10,10,5), c(0,0,1000,1000), col = blues[2], border = NA)
polygon(c(10,15,15,10), c(0,0,1000,1000), col = blues[3], border = NA)
polygon(c(15,20,20,15), c(0,0,1000,1000), col = blues[4], border = NA)
lines(days, 3*Hbs, lwd=2, col='red')

plot(1:length(drug_regimen)/24, out$retic_percent, type='l',ylim = c(0,10))
lines(days, Retics, lwd=2, col='red')

