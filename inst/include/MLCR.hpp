
#ifndef MLCR_hpp
#define MLCR_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type MLCR(objective_function<Type> *obj) {
  
  DATA_SCALAR(Linf);
  DATA_SCALAR(K);
  DATA_SCALAR(Lc);
  DATA_SCALAR(b);
  DATA_INTEGER(nbreaks);
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss); // Vector of sample sizes for mean lengths
  DATA_VECTOR(CPUE);
  DATA_INTEGER(CPUEisnormal); // 0 for lognormal, 1 for normal
  DATA_INTEGER(isWPUE); // 0 for NPUE, 1 for WPUE
  
  PARAMETER_VECTOR(Z);
  PARAMETER_VECTOR(yearZ);
  
  int count = Lbar.size();
  
  // Calculate a, v, r, s, denom, num, numsum, Lpred
  modelOutput<Type> output;
  if(nbreaks == 0) {
    output = model_output_eq(Z(0), Linf, K, Lc, b, count);
  } else {
    output = model_output(Z, yearZ, Linf, K, Lc, b, nbreaks, count);
  }
  
  vector<Type> Lpred(count);
  Lpred = output.Lpred;
  REPORT(Lpred);
    
  // Get population abundance (NPUE) or biomass (WPUE)
  vector<Type> population(count);
  vector<Type> Ipred;
  if(isWPUE == 0) {
    population = output.denom;
  } else {
    population = output.biomass;
  }
  
  // Analytical solution for q
  Type q = estimate_q(CPUE, population, CPUEisnormal, count);
  ADREPORT(q);
  
  // Calculate Ipred = qN or qB
  Ipred = calculate_Ipred(population, q, count);
  REPORT(Ipred);
  
  // Analytical solution for sigmaL
  Type sigmaL = estimate_sigmaL(Lbar, Lpred, ss, count);  
  ADREPORT(sigmaL);
  
  // Analytical solution for sigmaI
  Type sigmaI = estimate_sigmaI(CPUE, Ipred, count, CPUEisnormal);
  ADREPORT(sigmaI);
  
  // Initialize negative log-likelihood
  vector<Type> nllc(2);
  nllc.setZero();
  
  // Negative log-likelihood of mean lengths
  nllc(0) = nll_Lbar(Lbar, Lpred, ss, sigmaL, count);
  
  // Negative log-likelihood of catch rates
  nllc(1) = nll_CPUE(CPUE, Ipred, sigmaI, count, CPUEisnormal);
  
  Type nll = nllc.sum();

  REPORT(Z);
  REPORT(yearZ);
  REPORT(nllc);
  REPORT(nll);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

