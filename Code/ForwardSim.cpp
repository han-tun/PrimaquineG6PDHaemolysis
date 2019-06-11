#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// This function computes the effective dose given a dosing schedule vector
double compute_effective_dose(NumericVector drug_regimen, int t, int T_m){
  // drug_regimen: a vector of hourly primaquine `equivalent' doses 
  // (in practice only dosed once a day, this can later changed to drug concentration)
  // t: current time point
  // T_m: time delay to reach full effect
  if(t > drug_regimen.size()) throw std::range_error("Inadmissible value for t");
  double effective_dose = 0.0;
  if(t > 0){
    // index at which we start counting for the effective dose
    int jmin = max(IntegerVector::create(0,t-T_m)); 
    // sum over previous doses 
    for(int j=jmin; j<t; ++j) effective_dose += drug_regimen[j];
  }
  effective_dose = effective_dose / ((double) T_m); //normalise to get effective dose
  return effective_dose;
}

// [[Rcpp::export]]
// Convert an Hb value to a number of RBCs in the body
double Hb_to_NumberCells(double Hb, double BloodVolume, double MeanCellHb){
  double Ncells = pow(10,6) * Hb * BloodVolume * 10 / MeanCellHb;
  return Ncells;
}

// [[Rcpp::export]]
// Convert a number of RBCs in the body to an Hb concentration
double NumberCells_to_Hb(double Ct, double BloodVolume, double MeanCellHb){
  double Hb = pow(10,-6) * Ct * MeanCellHb / (10 * BloodVolume);
  return Hb;
}

// [[Rcpp::export]]
// Compute the fold change from baseline in the number of new normoblasts produced 
// given the number of cells in the body
double Fold_Change_Production(double rho_max,       // Max fold increase
                              double Cstar,         // Steady state number of cells
                              double C_t_minus_1,   // Total number of cells at previous timepoint
                              double C_50_rho       // mid point
){
  double rho;
  rho = rho_max/(1 + exp(log(rho_max-1)/(Cstar - C_50_rho)*(C_t_minus_1 - C_50_rho)));
  return rho;
}

// [[Rcpp::export]]
// Compute the steady state life span for a given effective dose
double compute_steady_state_life_span(double T_E_star,        // Lifespan when no drug
                                      double E_max,           // Maximum reduction in lifespan when given PMQ
                                      double effectivePMQdose,// current effective dose
                                      double x_50,            // mid-point
                                      double PMQ_slope){      // slope parameter
  double SteadyState_LifeSpan;
  SteadyState_LifeSpan = T_E_star*(1 - E_max/(1 + exp(-PMQ_slope * (log(effectivePMQdose) - log(x_50) ))) );
  return SteadyState_LifeSpan;
}

// [[Rcpp::export]]
// Compute the age at which retics enter the circulation
double compute_transit_time(double C_t_minus_1, // number of cells at previous time point
                            double Hb_50_circ,  // mid point (hb units)
                            double k,           // slope
                            double BloodVolume, // Volume of blood
                            double MeanCellHb){ // Hb per cell
  double transit;
  // next line should be changed: no need to recalculate each time
  double C_50_circ = Hb_to_NumberCells(Hb_50_circ,BloodVolume,MeanCellHb);
  transit = floor( (24*2.56)/(1+exp(-k * (C_t_minus_1 - C_50_circ) ))) + 24;
  return transit;
}

// [[Rcpp::export]]
// Count the number of retics in circulation
double CountRetics(int transit, NumericVector reticulocytes, NumericVector erythrocytes,
                   int extra_retics){
  double Total_retics = 0;
  int T_retic = reticulocytes.size();
  for(int i=transit; i<T_retic; ++i) Total_retics += reticulocytes[i];
  for(int i=0; i<extra_retics; ++i) Total_retics += erythrocytes[i];
  return Total_retics;
}

// [[Rcpp::export]]
// The runs a deterministic simulation
List forward_sim(NumericVector drug_regimen,  // dose present at each hour (mg/kg)
                 double rho_max,              // maximum fold increase in RBC production
                 double Hb_steady_state,      // steady state Hb (g/dL)
                 double Hb_50_rho,            // Hb corresponding to half of max fold increase (g/dL)
                 double Hb_50_circ,           // Hb corresponding to half transit time value (g/dL)
                 double k,                    // steepness of transit curve
                 double E_max,                // maximum relative decrease in length of RBC life span
                 double PMQ_slope,            // slope parameter of dose response curve
                 double x_50,                 // dose giving half decrease in life span (mg/kg)
                 int T_E_star,                // steady state duration of RBC life (hours)
                 int T_m,                     // delay in reaching effective dose (hours)
                 double MeanCellHb,           // mean corpuscular Hb (picogram)
                 double BloodVolume           // Total blood volume (liters)
){

  
  // Hack variable
  int extra_retics = 12; 
  // set up variables we need for the simulation
  int nhours = drug_regimen.size(); // This determines the length of the simulation in hours
  int i, t;                         // i iterates over age distributions, t over time of simulation
  int T_nmblast  = 24*5;            // number of hours of maturation of normoblasts
  int T_retic = 24*5;               // lifespan of a reticulocyte, including marrow and circulating periods
  int T_RBC_max = 24*150;           // max duration of RBC life: arbitrary to some extent (hours)
  double SS_LifeSpan;        // The drug dependent life span of an RBC (true continuous value)
  int T_lifeSpan;                         // The drug dependent life span of an RBC (rounded up to hours)
  double rho;                             // stores the fold change increase in production of normoblasts
  double Total_Eryths, Total_retics;      // Total circulating erythrocytes and retics
  double exponent_factor;                 // Used to compute the effect
  
  NumericVector effectiveDose(nhours);    // this is used for the lag effect, we can reduce down to single value later (currently tracking it)
  NumericVector Hb(nhours);               // The Hb at each point in the simulation
  NumericVector C_t(nhours);              // The number of circulating cells at each point in the simulation
  NumericVector retic_percent(nhours);    // The retic count (%) at each point in the simulation
  
  // vectors to store red blood cell age distributions
  NumericVector erythrocytes(T_RBC_max);  // Distribution of erythrocytes in circulation
  NumericVector normoblasts(T_nmblast);   // Distribution of normoblasts 
  NumericVector reticulocytes(T_retic);   // Distribution of reticulocytes 
  NumericVector temp_circ(T_RBC_max);     // for storing: temporary vectors
  NumericVector temp_normoblasts(T_nmblast);
  NumericVector temp_retics(T_retic);
  
  
  //***** some checks on input values ****
  if(min(drug_regimen) < 0) throw std::range_error("Inadmissible value");
  if(nhours == 0) return List::create(NA_INTEGER);
  if(T_E_star > T_RBC_max) throw std::range_error("T_E_star cannot be greater than T_RBC_max");
  // More thorough input checks are needed....
  
  // ***** Calculate initial values *****
  
  // The number of red blood cells in the body at steady state
  double Cstar = Hb_to_NumberCells(Hb_steady_state,BloodVolume,MeanCellHb);
  double C_50_rho = Hb_to_NumberCells(Hb_50_rho,BloodVolume,MeanCellHb);
  C_t[0] = Cstar; 
  effectiveDose[0] = 0.0;
  // The steady state lifespan: assuming at no drug steady state
  SS_LifeSpan = compute_steady_state_life_span(T_E_star,
                                               E_max,
                                               effectiveDose[0],
                                                            x_50,
                                                            PMQ_slope);
  T_lifeSpan = round(SS_LifeSpan); // the rounded up approximation
  Hb[0] = Hb_steady_state;
  
  // compute time at which retics are in the circulation
  int transit = compute_transit_time(C_t[0],
                                     Hb_50_circ,
                                     k, BloodVolume, MeanCellHb);
  // baseline production of normoblasts per hour
  double lambda = Cstar/( (double) T_lifeSpan + 24 + T_retic - transit);    
  
  // kill off 50% each day after T_E_star
  // initialise our parameters for distribution of RBCs in marrow and circulation
  for(i=0; i<T_RBC_max; i++) { erythrocytes[i] = lambda; }
  for(i = T_lifeSpan+1; i<T_RBC_max; ++i){
    exponent_factor = (double) (i - T_lifeSpan)/24.0;
    erythrocytes[i] *= pow(0.5,exponent_factor);
  }
  for(i=0; i<T_nmblast; i++) { normoblasts[i] = lambda; }
  for(i=0; i<T_retic; i++) { reticulocytes[i] = lambda; }
  
  Total_retics = CountRetics(transit, reticulocytes, erythrocytes,extra_retics);
  Total_Eryths = sum(erythrocytes);
  
  retic_percent[0] = 100 * Total_retics/(Total_retics + Total_Eryths);
  
  for(t=1; t < nhours; ++t){
    
    // Compute the multiplication factor of basal normoblast production
    rho = Fold_Change_Production(rho_max, Cstar, C_t[t-1], C_50_rho);
    
    // We move the RBCs from one compartment to the next
    temp_normoblasts = clone(normoblasts);
    normoblasts[0] = rho*lambda;      // the number of new normoblasts made at time t
    for(i = 1; i < T_nmblast; ++i){   // move all the normoblasts along by one
      normoblasts[i] = temp_normoblasts[i-1];
    }
    
    temp_retics = clone(reticulocytes);
    reticulocytes[0] = temp_normoblasts[T_nmblast-1];
    for(i = 1; i < T_retic; ++i){ // move all the retics along by one
      reticulocytes[i] = temp_retics[i-1];
    }
    
    // This part models the drug dependent killing as a shift in lifespan
    // Compute the life span at steady state for the current dose
    effectiveDose[t] = compute_effective_dose(drug_regimen, t, T_m);
    SS_LifeSpan = compute_steady_state_life_span(T_E_star,
                                                 E_max,
                                                 effectiveDose[t],
                                                              x_50,
                                                              PMQ_slope);
    T_lifeSpan = round(SS_LifeSpan);
    
    // update the age distribution of erythrocytes
    temp_circ[0] = erythrocytes[0];
    erythrocytes[0] = temp_retics[T_retic-1];
    for(i=1; i < T_RBC_max; ++i){
      temp_circ[i] = erythrocytes[i];
      if(i > T_lifeSpan) {
        // from T_lifeSpan onwards, reduce by half each day
        exponent_factor = (double) (i - T_lifeSpan)/24.0;
        erythrocytes[i] = pow(0.5, exponent_factor) * temp_circ[i-1];
      } else {
        erythrocytes[i] = temp_circ[i-1]; // move along by one
      }
    }
    
    // Count the number of retics and erythrocytes in circulation
    Total_retics = CountRetics(transit, reticulocytes, erythrocytes,extra_retics);
    Total_Eryths = sum(erythrocytes);
    
    // these values are stored in order to be returned at the end
    C_t[t] = Total_retics + Total_Eryths;
    Hb[t] = NumberCells_to_Hb(C_t[t],BloodVolume,MeanCellHb);
    retic_percent[t] = 100* Total_retics/(Total_retics + Total_Eryths);
    
    // calculate the updated Hb dependent transit time for the reticulocytes
    transit = compute_transit_time(C_t[t],
                                   Hb_50_circ,
                                   k, BloodVolume, MeanCellHb);
    
    
  }
  return List::create(Named("Hb", Hb),
                      Named("retic_percent", retic_percent),
                      Named("normoblasts", normoblasts),
                      Named("reticulocytes", reticulocytes),
                      Named("erythrocytes", erythrocytes),
                      Named("SS_LifeSpan", SS_LifeSpan),
                      Named("Retic_Transit", transit),
                      Named("T_retic", T_retic),
                      Named("effectiveDose", effectiveDose));
  
}