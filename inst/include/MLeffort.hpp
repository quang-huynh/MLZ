
#ifndef MLeffort_hpp
#define MLeffort_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type MLeffort(objective_function<Type> *obj) {
  
  DATA_SCALAR(Linf);
  DATA_SCALAR(K);
  DATA_SCALAR(a0);
  DATA_SCALAR(Lc);
  
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss);
  DATA_VECTOR(eff);
  
  DATA_SCALAR(eff_init);
  DATA_INTEGER(n_age);
  DATA_INTEGER(n_season);
  DATA_INTEGER(obs_season);
  
  DATA_SCALAR(timing);
  DATA_INTEGER(logpar);
  
  PARAMETER(logq);
  PARAMETER(logM);
    
  int n_yr = Lbar.size();
  int y;
  int a;
  int obs_ind = obs_season - 1;
  int astep = n_age * n_season;
  
  Type q = 0.;
  Type M = 0.;
  
  if(logpar == 0) {
    q = logq;
    M = logM;
  }
  if(logpar == 1) {
    q = exp(logq);
    M = exp(logM);
  }
  ADREPORT(q);
  ADREPORT(M);
  
  Type Z_init = q * eff_init + M;
  vector<Type> Z(n_yr);
  
  vector<matrix<Type> > N(n_yr);
  //array<Type> N(n_yr,astep,n_season);
  matrix<Type> Nobs(n_yr,astep);
  vector<Type> La(astep);
  vector<Type> age(astep);
  
  vector<Type> Lpred(n_yr);
  
  Type seasD = n_season;
  Type ac = a0 - log(1 - Lc/Linf)/K;
  for(a=0;a<astep;a++) {
    Type ageD = a;
    age(a) = ac + ageD/seasD;
    La(a) = Linf * (1 - exp(-K*(age(a) - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;
  REPORT(Z);
  
  N(0) = ML_effort_Neq(astep, n_season, Z_init, Z(0), seasD);
  
  for(y=1;y<n_yr;y++) {
    N(y) = ML_effort_N(astep, n_season, Z(y-1), Z(y), seasD, N(y-1));
  }
  
  for(y=0;y<n_yr;y++) {
    Type num = 0.;
    Type denom = 0.;
    for(a=0;a<astep;a++) {
      Nobs(y,a) = N(y)(a,obs_ind) * exp(-Z(y) * timing/seasD);
      num += Nobs(y,a) * La(a);
      denom += Nobs(y,a);
    }
    Lpred(y) = num/denom;
  }
  REPORT(Lpred);
  
  Type sigma = estimate_sigmaL(Lbar, Lpred, ss, n_yr);
  ADREPORT(sigma);
  
  Type nll = nll_Lbar(Lbar, Lpred, ss, sigma, n_yr); 
  
  return nll;
  
}
  
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  
#endif
