  
  real<lower=0> EPS;      // correction to avoid neg. values
  real<lower=0> abs_tol;  // absolute tolerance for BDF
  real<lower=0> rel_tol;  // relative tolerance for BDF
  int<lower=1> max_iter;  // max number of iterations for BDF
  
  
  int K;                   // number of age classes
  vector[K] age_dist;
  int pop_t;               // total population
  real tswitch;
  
  // Data to fit
  int D;                   // number of days with reported incidence
  int incidence_cases[D];  // overal incidence for W weeks
  int incidence_deaths[D]; // overal incidence for W weeks
  int agedistr_cases[K];   // number of cases at tmax for the K age classes
  int agedistr_deaths[K];  // mortality at tmax for the K age classes
  
  // Prior hyperparameters
  real p_beta;
  real p_eta[2];
  real p_pi[2];
  real p_epsilon[2];
  real p_rho[2];
  real p_phi;
  real p_xi;
  real p_nu;
  real p_psi[2];
  
  // Fixed parameters
  real p_incubation;
  real p_infectious;
  int G;
  real p_gamma[G];
  
  // Simulation
  real t0;           // starting time
  int t_data;        // time of first data
  int S;
  real ts[S];        // time bins
  real contact[K*K];
