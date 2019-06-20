#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// This function computes the effective dose given a dosing schedule vector
List compute_effective_dose(NumericVector drug_regimen, int t, double lambda_memory = -2){
  // drug_regimen: a vector of hourly primaquine `equivalent' doses 
  // (in practice only dosed once a day, this can later changed to drug concentration)
  // t: current time point
  if(t > drug_regimen.size()) throw std::range_error("Inadmissible value for t");
  double effective_dose = 0.0;
  NumericVector weights(t); 
  double dt = 1.0 / 24.0;
  double proportional_sum;
  if(t > 0){
    proportional_sum = 1.0 - exp(lambda_memory*t/24);
    for(int i=0; i < t; ++i) weights[i] = dt * exp(lambda_memory * (double) (t-i)/24 );
    weights = proportional_sum * weights/ sum(weights);
    effective_dose = sum(weights * drug_regimen);
  }
  //return effective_dose;
  return List::create(Named("weights", weights),
                      Named("pp", proportional_sum),
                      Named("effectiveDose", effective_dose));
}