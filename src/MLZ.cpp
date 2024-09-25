#define TMB_LIB_INIT R_init_MLZ
#include <TMB.hpp>
#include <tiny_ad/beta/pbeta.hpp>
#include "../inst/include/functions.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  
  if(model == "ML") {
    return ML(this);
  } else if(model == "MLCR") {
    return MLCR(this);
  } else if(model == "MSM1S") {
    return MSM1S(this); 
  } else if(model == "MSM23") {
    return MSM23(this);
  } else if(model == "MLeffort") {
    return MLeffort(this);
  }
  return 0;
}
