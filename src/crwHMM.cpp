
#include <TMB.hpp>
#include <cmath>

using namespace density;
using std::sqrt;


// Function for creating cov mat for ctcrw SSM

// Function for creating projection matrix for ctcrw SSM

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA
  DATA_MATRIX(Y);	            //  (x, y) observations
  DATA_VECTOR(dt);         	//  time diff in some appropriate unit. this should contain dt for both interp and obs positions.
  DATA_VECTOR(mu0);        //  initial state mean
  DATA_SCALAR(V0);          // Prior variance of first location
  DATA_IVECTOR(isd);          //  indexes observations vs. interpolation points
  DATA_IVECTOR(obs_mod); //  indicates which obs error model to be used
  
  DATA_INTEGER(nbStates);
  
  // for KF observation model
  DATA_VECTOR(m);             //  m is the semi-minor axis length
  DATA_VECTOR(M);             //  M is the semi-major axis length
  DATA_VECTOR(c);             //  c is the orientation of the error ellipse
  // for LS observation model
  DATA_MATRIX(K);                 // error weighting factors for LS obs model
  
  // PROCESS PARAMETERS
  // for CRW
  PARAMETER_VECTOR(log_beta);				
  PARAMETER_VECTOR(log_sigma);
  // random variables
  PARAMETER_MATRIX(mu); /* State matrix */
  
  // OBSERVATION PARAMETERS
  // for KF OBS MODEL
  PARAMETER(l_psi); 				// error SD scaling parameter to account for possible uncertainty in Argos error ellipse variables
  // for LS OBS MODEL
  PARAMETER_VECTOR(l_tau);     	// error dispersion for LS obs model (log scale)
  PARAMETER(l_rho_o);             // error correlation
  
  // Transform parameters
  vector<Type> beta = exp(log_beta);
  vector<Type> sigma = exp(log_sigma);
  Type psi = exp(l_psi);
  vector<Type> tau = exp(l_tau);
  Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  
  int timeSteps = dt.size();
  
  /* Define likelihood */
  //parallel_accumulator<Type> jnll(this);
  matrix<Type> mu_bar_i(2,nbStates);
  vector<Type> sig_i;
  Type nll_proc=0.0;
  Type nll_obs=0.0;
  Type nll=0.0;
  vector<Type> x;
  
  // CRW
  // Setup object for evaluating multivariate normal likelihood
  matrix<Type> cov_obs(2, 2);
  cov_obs.setZero();
  
  //initial locs
  // time 0
  nll_proc -= dnorm(mu(0,0), mu0(0), V0,  1) + 
    dnorm(mu(1,0), mu0(1), V0,  1); 
  //time 1
  for(int j=0; j<nbStates; j++) {
    sig_i = dt(1)*sigma/(2*beta);
  }
  nll_proc -= dnorm(mu(0,1), mu(0,0), sig_i,  1) + 
    dnorm(mu(1,1), mu(1,0), sig_i,  1); 

  for(int i = 2; i < timeSteps; i++) {
    for(int j=0; j<nbStates; j++) {
      mu_bar_i.col(j) = (1+dt(i)/dt(i-1)-beta(j)*dt(i))*mu.col(i-1) + 
        (beta(j)*dt(i)-dt(i)/dt(i-1))*mu.col(i-2);
      sig_i(j) = dt(i)*sqrt(dt(i-1))*sigma(j);
    }
    nll_proc -= dnorm(mu(0,i), mu_bar_i(0), sig_i, 1) + 
      dnorm(mu(1,i), mu_bar_i(1), sig_i, 1);
  }
  
  // OBSERVATION MODEL
  // 2 x 2 covariance matrix for observations
  
  for(int i = 0; i < timeSteps; ++i) {
  if(isd(i) == 1) {
    if(obs_mod(i) == 0) {
      // Argos Least Squares observations
      Type s = tau(0) * K(i,0);
      Type q = tau(1) * K(i,1);
      cov_obs(0,0) = s * s;
      cov_obs(1,1) = q * q;
      //cov_obs(0,1) = s * q * rho_o;
      //cov_obs(1,0) = cov_obs(0,1);
    } else if(obs_mod(i) == 1) {
      // Argos Kalman Filter (or Kalman Smoothed) observations
      Type z = sqrt(Type(2.0));
      Type s2c = sin(c(i)) * sin(c(i));
      Type c2c = cos(c(i)) * cos(c(i));
      Type M2  = (M(i) / z) * (M(i) / z);
      Type m2 = (m(i) * psi / z) * (m(i) * psi / z);
      cov_obs(0,0) = (M2 * s2c + m2 * c2c);
      cov_obs(1,1) = (M2 * c2c + m2 * s2c);
      cov_obs(0,1) = (0.5 * (M(i) * M(i) - (m(i) * psi * m(i) * psi))) * cos(c(i)) * sin(c(i));
      cov_obs(1,0) = cov_obs(0,1);
    }
    nll_obs +=  MVNORM(cov_obs)(Y.col(i) - mu.col(i));   // CRW innovations


  }
  }
  
  nll = nll_proc + nll_obs;
  
  ADREPORT(beta);
  ADREPORT(sigma);
  ADREPORT(rho_o);
  ADREPORT(tau);
  ADREPORT(psi);
  
  REPORT(nll_proc);
  REPORT(nll_obs);
  
  return nll;
}

