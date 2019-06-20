library(truncnorm)

estimate_sim_profile = function(params,drug_regimen){
  out = forward_sim(drug_regimen = as.double(drug_regimen), # mg/kg doses
                    rho_max = params['rho_max'],
                    Hb_steady_state = params['Hb_star'],
                    Hb_50_rho = params['Hb_50_rho'],
                    Hb_50_circ = params['Hb_50_circ'],
                    k = 10^(-6),
                    E_max = params['E_max'],
                    x_50 = params['x_50'],
                    T_E_star = round(params['T_E_star']),
                    PMQ_slope = params['PMQ_slope'],
                    MeanCellHb = 30,
                    BloodVolume = 5)
  return(out)
}

plot_Hb_data = function(Haemodata_Analysis, col_line=NA){
  mean_Hbs = array(dim=length(unique(Haemodata_Analysis$Day)))
  for(dd in unique(Haemodata_Analysis$Day)){
    ind = Haemodata_Analysis$Day==dd
    points(dd-1, mean(Haemodata_Analysis$Hb[ind]), pch=18, col='red')
    if(is.na(col_line)) col_line='black'
    lines(c(dd,dd)-1, quantile(Haemodata_Analysis$Hb[ind], probs = c(.1,.9)),col=col_line)
    mean_Hbs[which(dd==unique(Haemodata_Analysis$Day))] = mean(Haemodata_Analysis$Hb[ind])
  }
  lines(unique(Haemodata_Analysis$Day)-1, mean_Hbs, col='red', lwd=2)
}

plot_Retic_data = function(Haemodata_Analysis, col_line=NA){
  mean_Retics = array(dim=length(unique(Haemodata_Analysis$Day)))
  for(dd in unique(Haemodata_Analysis$Day)){
    ind = Haemodata_Analysis$Day==dd
    if(is.na(col_line)) col_line='black'
    points(dd-1, mean(Haemodata_Analysis$Rectic_percent[ind]), pch=18, col='red')
    lines(c(dd,dd)-1, quantile(Haemodata_Analysis$Rectic_percent[ind], probs = c(.1,.9),na.rm = T),col=col_line)
    mean_Retics[which(dd==unique(Haemodata_Analysis$Day))] = mean(Haemodata_Analysis$Rectic_percent[ind])
  }
  lines(unique(Haemodata_Analysis$Day)-1, mean_Retics, col='red', lwd=2)
}

generate_from_prior = function(n=1){
  params_list = list()
  Hb_star = 15
  for(i in 1:n){
    params = list(rho_max = rtruncnorm(n = 1, mean = 5, sd = 1, a = 2),
                  Hb_50_rho = rtruncnorm(n = 1, mean = 10, sd = 2, a = 0, b = Hb_star),
                  Hb_50_circ = rtruncnorm(n = 1, mean = 10, sd = 2, a = 0, b = Hb_star),
                  E_max = rbeta(n = 1, shape1 = 10, shape2 = 20),
                  x_50 = rtruncnorm(n = 1, mean = 15, sd = 10, a = 5)/60,
                  T_E_star = rtruncnorm(n = 1, mean = 100, sd = 30, a=50, b=120)*24,
                  PMQ_slope = rexp(n = 1,rate = 1),
                  Hb_star = Hb_star)
    params = unlist(params)
    params_list[[i]] = params
  }
  
  return(params_list)
}

prior_density = function(params){
  pi_p = log(dtruncnorm(x = params['rho_max'], mean = 6, sd = 1,a=2)) +
    log(dtruncnorm(x = params['Hb_50_rho'],mean = 10, sd = 5, a = 0, b = params['Hb_star'])) +
    log(dtruncnorm(x = params['Hb_50_circ'],mean = 10, sd = 5, a = 0, b = params['Hb_star'])) +
    dbeta(x = params['E_max'], shape1 = 10, shape2 = 20,log = T) +
    log(dtruncnorm(x = params['x_50']*60, mean = 15, sd = 10, a = 5)) + 
    log(dtruncnorm(x = params['T_E_star']/24, mean = 100, sd = 10, a = 50,b=120)) +
    dexp(x = params['PMQ_slope'], rate = 1,log = T)
  
  return(pi_p)
}




compute_post_density = function(params, drug_regimen,Haemodata_Analysis){
  pi_p = prior_density(params)
  if(!is.infinite(pi_p)){
    out = forward_sim(drug_regimen = as.double(drug_regimen), # mg/kg doses
                      rho_max = params['rho_max'],
                      Hb_steady_state = params['Hb_star'],
                      Hb_50_rho = params['Hb_50_rho'],
                      Hb_50_circ = params['Hb_50_circ'],
                      E_max = params['E_max'],
                      x_50 = params['x_50'],
                      T_E_star = round(params['T_E_star']),
                      PMQ_slope = params['PMQ_slope'])
    ind_out_days = which((0:(length(out$Hb)-1) %% 24) == 0)
    log_post = pi_p + 
      sum(dnorm(x = out$Hb[ind_out_days[Haemodata_Analysis$Day]] - Haemodata_Analysis$Hb, mean = 0, sd = 1, log = T)) +
      sum(dcauchy(x = out$retic_percent[ind_out_days[Haemodata_Analysis$Day]] - Haemodata_Analysis$Rectic_percent,
                  location = 0, scale = 1.5, log = T), na.rm = T)
    return(log_post)
  } else {
    return(-Inf)
  }
}



plot_dose_response = function(params, col = 'black',lwd=1,add=F,ylim=NA){
  ds = seq(0,45,length.out = 100)/60
  ssls = array(dim = length(ds))
  for(i in 1:length(ds)){
    ssls[i] = compute_steady_state_life_span(T_E_star = params['T_E_star'],
                                             E_max = params['E_max'], # fraction decrease
                                             effectivePMQdose = ds[i],
                                             x_50 = params['x_50'],
                                             PMQ_slope = params['PMQ_slope'])
  }
  if(add){
    lines(ds*60, ssls/24,lwd=lwd,col=col)
  } else {
    if(is.na(ylim)) ylim=range(ssls/24)
    plot(ds*60, ssls/24, type='l',xlab = 'Effective PMQ dose in 60 kg adult',
         ylab = 'Days survival of erythrocytes',lwd=lwd,col=col,ylim=ylim)
    title('Dose-response curve')
  }
  
}

