#define TMB_LIB_INIT R_init_MLZ
#include <TMB.hpp>
#include <tiny_ad/beta/pbeta.hpp>
#include <cppad/cppad.hpp>
#include "../inst/include/functions.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);
  
  if(model == "ML") {
      #include "../inst/include/ML.h"
  } else 
      if(model == "MLCR") {
          #include "../inst/include/MLCR.h" 
      } else
          if(model == "MSM1S") {
              #include "../inst/include/MSM1S.h" 
          } else 
              if(model == "MSM23") {
                  #include "../inst/include/MSM23.h" 
              } else
                  if(model == "MLeffort") {
                      #include "../inst/include/MLeffort.h" 
                  } else
                      error("No model found.");
  
  return 0;
}
