#include <Rcpp.h>
using namespace Rcpp;

double compute_effect(double time_at_dose,
                      double time_to_max_effect){
  double effect;
  // calculate the effect based on the time at dose and the previous dose
  effect = pow(min(NumericVector::create(time_at_dose/time_to_max_effect,1)), 2);
  return effect;
}

double compute_steady_state_life_span(double LifeSpan_norm,
                                      double LifeSpan_min,
                                      double PMQdose,
                                      double Effect50_dose,
                                      double PMQ_slope){
  
  double SteadyState_LifeSpan;
  SteadyState_LifeSpan = (double) LifeSpan_norm - (LifeSpan_norm - LifeSpan_min)/(1 + exp(-PMQ_slope * (log(PMQdose) - log(Effect50_dose))));
  return SteadyState_LifeSpan;
}

double compute_transit_time(double HCT,
                            double mid_transit,
                            double hill_coef){
  double transit;
  transit = round( (24*2.56)/(1+exp( (mid_transit - HCT)/ hill_coef ))) + 24;
  return transit;
}

// [[Rcpp::export]]
List forward_sim(NumericVector drug_regimen,  // the mg/kg dose present at each hour
                 double max_production,       // maximum fold increase in RBC production
                 double mid_haematocrit,      // HCT corresponding to half of max fold increase
                 double HCT_steady_state,     // steady state HCT
                 double mid_transit,          // HCT corresponding to half transit time value
                 double hill_coef,            // steepness of transit curve
                 double time_to_max_effect,             // time lag value
                 double LifeSpan_min,            // Min length of RBC life span under max doses of PMQ
                 double PMQ_slope,            // slope parameter of dose response curve
                 double Effect50_dose,        // Dose giving half decrease in life span in dosing units
                 int LifeSpan_norm          // normal duration of RBC life (hours)
){
  
  
  int nhours = drug_regimen.size(); // This determines the length of the simulation
  // some checks on input values
  if(min(drug_regimen) < 0) throw std::range_error("Inadmissible value");
  if(nhours == 0) return List::create(NA_INTEGER);
  
  int i, t;                         // i iterates over age distributions, t over time of simulation
  int T_nmblast  = 24*5;            // number of hours of maturation of normoblasts
  int T_retic = 24*5;               // lifespan of a reticulocyte, including marrow and circulating periods
  int T_RBC_max = 24*125;           // max duration of RBC life (hours)
  // check range of LifeSpan_norm
  if(LifeSpan_norm > T_RBC_max) throw std::range_error("LifeSpan_norm cannot be greater than T_RBC_max");
  
  int T_lifeSpan;                   // The drug dependent life span of an RBC
  NumericVector SS_LifeSpan(nhours); 
  double Current_LifeSpan;
  
  double lambda  = 1e5;             // Arbitrary baseline production of cells/kg/hour
  
  // slope parameter for the normoblast production function, dependent on input arguments
  double slope_k = log(max_production -  1)/(HCT_steady_state - mid_haematocrit);    
  
  double rho = 1;                   // stores the fold change increase in production of normoblasts
  double frac_change = 1.0;

  double Total_Eryths, Total_retics;      // Total circulating erythrocytes and retics
  NumericVector haematocrit(nhours);      // The haematocrit at each point in the simulation
  NumericVector retic_percent(nhours);    // The retic count (%) at each point in the simulation
  
  // Red blood cell distributions
  NumericVector erythrocytes(T_RBC_max);      // Distribution of erythrocytes in circulation
  NumericVector normoblasts(T_nmblast);       // Distribution of normoblasts 
  NumericVector reticulocytes(T_retic);       // Distribution of reticulocytes 
  NumericVector temp_circ(T_RBC_max);         // for storing: temporary vectors
  NumericVector temp_marrow(T_nmblast);
  NumericVector temp_retics(T_retic);
  
  double exponent_factor;                 // Used to compute the effect
 
  int time_at_dose = 0;
  NumericVector effect(nhours);     // this is used for the lag effect
  effect[0] = compute_effect((double) time_at_dose, (double) time_to_max_effect);
  SS_LifeSpan[0] = compute_steady_state_life_span(LifeSpan_norm,
                                                  LifeSpan_min,
                                                  0.0,
                                                  Effect50_dose,
                                                  PMQ_slope);
  T_lifeSpan = round(SS_LifeSpan[0]);
  // kill off 50% each day after LifeSpan_norm
  // initialise our parameters for distribution of RBCs in marrow and circulation
  // Introduce some random variation to smooth output
  erythrocytes = lambda * rnorm(T_RBC_max, 1, 0.01);
  for(i = T_lifeSpan+1; i<T_RBC_max; ++i){
    exponent_factor = (double) (i - T_lifeSpan)/24.0;
    erythrocytes[i] = pow(0.5,exponent_factor) * erythrocytes[i];
  }
  normoblasts = lambda * rnorm(T_nmblast, 1, 0.01);
  reticulocytes = lambda * rnorm(T_retic, 1, 0.01);
  
  // compute time at which retics are in the circulation
  int transit = compute_transit_time(haematocrit[0],
                                     mid_transit,
                                     hill_coef);
  
  Total_retics = 0;
  for(i=transit; i<T_retic; ++i) Total_retics += reticulocytes[i];
  Total_Eryths = sum(erythrocytes);
  double N_RBC_init = Total_Eryths + Total_retics; // number of circulating RBCs at steady state
  
  // this values are stored in order to be returned at the end
  haematocrit[0] = HCT_steady_state;
  retic_percent[0] = 100 * Total_retics/(Total_retics + Total_Eryths);
  
  double Previous_Dose = 0.0;         
  
  // Simulation up to nhours dependent on the primaquine dosing defined by drug_regimen[t]
  for(t=1; t < nhours; ++t){
    
    // Compute the multiplication factor of basal normoblast production
    rho = max_production/(1 + exp(slope_k*(haematocrit[t-1] - mid_haematocrit)));
    
    temp_marrow[0] = normoblasts[0];
    normoblasts[0] = rho*lambda;      // the number of new normoblasts made at time t
    for(i = 1; i < T_nmblast; ++i){   // move all the normoblasts along by one
      temp_marrow[i] = normoblasts[i];
      normoblasts[i] = temp_marrow[i-1];
    }
    
    temp_retics[0] = reticulocytes[0];
    reticulocytes[0] = temp_marrow[T_nmblast-1];
    for(i = 1; i < T_retic; ++i){ // move all the retics along by one
      temp_retics[i] = reticulocytes[i];
      reticulocytes[i] = temp_retics[i-1];
    }
    
    // This part models the drug dependent killing as a shift in lifespan
    // Compute the life span at steady state for the current dose
    SS_LifeSpan = compute_steady_state_life_span(LifeSpan_norm,
                                                 LifeSpan_min,
                                                 drug_regimen[t-1],
                                                 Effect50_dose,
                                                 PMQ_slope);
    
    if(drug_regimen[t] > 0.0 & time_at_dose < (time_to_max_effect+1)) time_at_dose++;
    
    // First we model the lag time in life span decrease
    if(Previous_Dose < drug_regimen[t]){ // increase in dosing
      if(Previous_Dose > 0.0){
        frac_change = min(NumericVector::create((drug_regimen[t]-Previous_Dose)/Previous_Dose, 1));
      } else {
        frac_change = 1.0;
      }
      time_at_dose = round(min(NumericVector::create(time_at_dose/time_to_max_effect, 1-frac_change))*time_to_max_effect);
    } 
    if(drug_regimen[t] > 0.0) {
      effect[t] = pow(min(NumericVector::create(time_at_dose/time_to_max_effect,1)), 2);
    } else {
      effect[t] = 0.0;
    }
    T_lifeSpan = round(LifeSpan_norm - (LifeSpan_norm-SS_LifeSpan)*effect[t]);
    
    //update the age distribution of erythrocytes
    temp_circ[0] = erythrocytes[0];
    erythrocytes[0] = temp_retics[T_retic-1];
    for(i=1; i < T_RBC_max; ++i){
      temp_circ[i] = erythrocytes[i];
      if(i > T_lifeSpan) {
        // from T_lifeSpan onwards, reduce by half each day
        exponent_factor = (double) (i - T_lifeSpan)/24.0;
        erythrocytes[i] = pow(0.5, exponent_factor) * temp_circ[i-1];
      } else {
        erythrocytes[i] = temp_circ[i-1];
      }
    }
    
    Total_retics = 0;
    for(i=transit; i<T_retic; ++i) Total_retics += reticulocytes[i];
    Total_Eryths = sum(erythrocytes);
    
    // this values are stored in order to be returned at the end
    haematocrit[t] = HCT_steady_state*(Total_Eryths+Total_retics)/(N_RBC_init);
    retic_percent[t] = 100* Total_retics/(Total_retics + Total_Eryths);
    
    // calculate the updated haematocrit dependent transit time for the reticulocytes
    transit = round( (24*2.56)/(1+exp( (mid_transit - haematocrit[t])/ hill_coef ))) + 24;
    Previous_Dose = drug_regimen[t];  
  }
  
  // Store all the circulating cells for plotting after (CAN REMOVE IN LATER VERSION)
  NumericVector CirculatingRBCs(T_RBC_max+T_retic);
  for(i=0; i<transit; ++i) CirculatingRBCs[i] = 0;
  for(i = transit; i<T_retic; ++i) CirculatingRBCs[i] = reticulocytes[i];
  for(i = 0; i<T_RBC_max; ++i) CirculatingRBCs[i+T_retic] = erythrocytes[i];
  
  return List::create(Named("haemtocrit", haematocrit),
                      Named("retic_percent", retic_percent),
                      Named("normoblasts", normoblasts),
                      Named("reticulocytes", reticulocytes),
                      Named("erythrocytes", erythrocytes),
                      Named("CiculatingRBCs", CirculatingRBCs),
                      Named("TLIFESPAN", T_lifeSpan),
                      Named("Retic_Transit", transit),
                      Named("T_retic", T_retic),
                      Named("effect", effect),
                      Named("frac_change", frac_change));
}