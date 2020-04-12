
  real<lower=0,upper=1> beta;             // base transmission rate
  real<lower=0,upper=1> eta;              // reduct. in transm. rate after quarant. measures
  vector<lower=0,upper=1> [K] epsilon;    // age-dependent mortality probability
  vector<lower=0,upper=1> [K-1] raw_rho;  // age-dependent reporting probability
  real<lower=0, upper=1> pii;             // number of cases at t0
  real<lower=0> phi[2];                   // variance parameters
  real<lower=0,upper=1> xi_raw;           // slope of quarantine implementation
  real<lower=0> nu;                       // shift of quarantine implementation
  real<lower=0,upper=1> psi;              // proportion of symptomatics
