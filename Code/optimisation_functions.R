optimise_regimen = function(params, min_total_dose = 360, min_increment = 2.5, duration=20,max_single_dose=45){
  dosing_vector = rep(5,duration)
  total_dose = sum(dosing_vector)
  while(total_dose < min_total_dose){
    dosing_vector = optimal_add(min_increment, dosing_vector,params,max_single_dose)
    total_dose = sum(dosing_vector)
  }
  return(dosing_vector)
}

optimal_add = function(min_increment, dosing_vector,params, max_single_dose){
  deltas = array(NA, length(dosing_vector))
  for(j in length(dosing_vector):1){
    test_dosing_vector = dose_increment(dosing_vector,j,min_increment,max_single_dose)
    if(any(is.na(test_dosing_vector))) {
      deltas[j] = NA
    } else {
      drug_regimen = array(sapply(test_dosing_vector, rep, 24))/60
      out = estimate_sim_profile(params = params,drug_regimen = drug_regimen)
      deltas[j] = min(diff(out$Hb[seq(1,length(out$Hb),by=24)]))
    }
    
  }
  dosing_vector = dose_increment(dosing_vector,which.max(deltas),min_increment,max_single_dose)
  return(dosing_vector)
}

dose_increment = function(dosing_vector,j,min_increment,max_single_dose){
  # this makes sure the regimen is ascending
  if(j<1 | j>length(dosing_vector)) stop('wrong index j : dose_increment function')
  if(dosing_vector[j] >= max_single_dose) return(NA)
  dosing_vector[j] = dosing_vector[j]+min_increment
  dosing_vector[1:length(dosing_vector)>j & dosing_vector<dosing_vector[j]] = dosing_vector[j]
  return(dosing_vector)
}
