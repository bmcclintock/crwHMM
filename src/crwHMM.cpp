
#include <TMB.hpp>
#include <cmath>

using namespace density;
using std::sqrt;

#define LOGZERO NAN;

/* implement the vector - matrix product */
template<class Type>
vector<Type> multvecmat(vector<Type>  A, matrix<Type>  B) {
  int nrowb = B.rows();
  int ncolb = B.cols(); 
  vector<Type> C(ncolb);
  for (int i = 0; i < ncolb; i++)
  {
    C(i) = Type(0);
    for (int k = 0; k < nrowb; k++){
      C(i) += A(k)*B(k,i);
    }
  }
  return C;
}

/* vector - matrix multiplication using atomic::matmal */
template<class Type>
vector<Type> multiply(vector<Type> x, matrix<Type> y){
  return atomic::matmul(matrix<Type>(x.matrix()),y).vec();
}

template<class Type>
Type eexp(Type x) {
  Type C = 0.0;
  if(!isnan(x)) C = exp(x);
  return C;
}

template<class Type>
Type eln(Type x) {
  Type C;
  if(x>0) C = log(x);
  else C = LOGZERO;
  return C;
}

template<class Type>
Type elnsum(Type x, Type y) {
  Type C;
  if(isnan(x) | isnan(y)){
    if(isnan(x)) C = y;
    else C = x;
  } else {
    if(x>y) C = x + eln(Type(1.0)+exp(y-x));
    else C = y + eln(Type(1.0)+exp(x-y));
  }
  return C;
}

template<class Type>
Type elnproduct(Type x, Type y) {
  Type C;
  if(isnan(x) | isnan(y)){
    C = LOGZERO;
  } else {
    C = x + y;
  }
  return C;
}

template <class Type>
int which_max(vector<Type> x){
  int first = 0;
  int last = x.size()-1;
  if (first==last) return last;
  int largest = first;
  
  while (++first<last+1)
    if (x(largest)<x(first))    
      largest=first;
    return largest;
}


/* Numerically stable forward algorithm based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf */
template<class Type>
Type forward_alg(vector<Type> delta, matrix<Type> trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = trMat.cols();
  Type logalpha;
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  vector<Type> lalpha(nbStates);
  vector<Type> lnewalpha(nbStates);
  Type sumalpha = LOGZERO;
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = LOGZERO;
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
      ldeltaG(j) = elnsum(ldeltaG(j),elnproduct(eln(delta(i)),ltrMat(i,j)));
    }
    lalpha(j) = elnproduct(ldeltaG(j),lnProbs(0,j));
    sumalpha  = elnsum(sumalpha,lalpha(j));
  }
  Type jnll = -sumalpha;
  lalpha -= sumalpha;
  for(int t=1; t < nbSteps; t++){
    sumalpha = LOGZERO;
    for(int j=0; j < nbStates; j++){
      logalpha = LOGZERO;
      for(int i=0; i < nbStates; i++){
        logalpha = elnsum(logalpha,elnproduct(lalpha(i),ltrMat(i,j)));
      }
      lnewalpha(j) = elnproduct(logalpha,lnProbs(t,j));
      sumalpha  = elnsum(sumalpha,lnewalpha(j));
    }
    jnll -= sumalpha;
    lalpha = lnewalpha - sumalpha;
  }
  return jnll;
}

/* forward algorithm using only scaling; this blows up when calculating exp(lnProbs) */
template<class Type>
Type forward_alg_scale(vector<Type> delta, matrix<Type> trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = trMat.cols();

  Type sumalpha;
  Type jnll = Type(0.0);
  vector<Type> Probs(nbStates);
  vector<Type> alpha = delta;
  
  for(int t=0; t < nbSteps; t++){
    Probs = exp(vector<Type>(lnProbs.row(t)));
    alpha = multvecmat(alpha,trMat) * Probs;
    sumalpha = alpha.sum();
    jnll -= log(sumalpha);
    alpha /= sumalpha;
  }
  return jnll;
}

/* Numerically stable viterbi based on Rabiner 1989 IEEE 77(2):257-286 */
template<class Type>
vector<int> viterbi(vector<Type> delta, matrix<Type> trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = trMat.cols();
  matrix<Type> phi(nbSteps,nbStates);
  matrix<int> psi(nbSteps,nbStates);
  vector<Type> tmpphi(nbStates);
  vector<int> states(nbSteps);
  states.setOnes();
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = LOGZERO;
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
      ldeltaG(j) = elnsum(ldeltaG(j),elnproduct(eln(delta(i)),ltrMat(i,j)));
    }
    psi(0,j) = 0;
    phi(0,j) = elnsum(ldeltaG(j),lnProbs(0,j));
  }
  for(int t=1; t < nbSteps; t++){
    for(int j=0; j < nbStates; j++){
      for(int i=0; i < nbStates; i++){
        tmpphi(i) = phi(t-1,i) + ltrMat(i,j);
      }
      psi(t,j) = which_max(tmpphi);
      phi(t,j) = tmpphi(psi(t,j)) + lnProbs(t,j);
    }
  }
  states(nbSteps-1) += which_max(vector<Type>(phi.row(nbSteps-1)));
  for(int t=0; t<(nbSteps-1); t++){
    states(nbSteps-2-t) += psi(nbSteps-1-t,states(nbSteps-1-t)-1);
  }
  return states;
}

/* Numerically stable forward log-probabilities based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf */
template<class Type>
matrix<Type> logAlpha(vector<Type> delta, matrix<Type> trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = trMat.cols();
  Type logalpha;
  matrix<Type> elnalpha(nbSteps,nbStates);
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = LOGZERO;
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
      ldeltaG(j) = elnsum(ldeltaG(j),elnproduct(eln(delta(i)),ltrMat(i,j)));
    }
    elnalpha(0,j) = elnproduct(ldeltaG(j),lnProbs(0,j));
  }
  for(int t=1; t < nbSteps; t++){
    for(int j=0; j < nbStates; j++){
      logalpha = LOGZERO;
      for(int i=0; i < nbStates; i++){
        logalpha = elnsum(logalpha,elnproduct(elnalpha(t-1,i),ltrMat(i,j)));
      }
      elnalpha(t,j) = elnproduct(logalpha,lnProbs(t,j));
    }
  }
  return elnalpha;
}

/* Numerically stable backward log-probabilities based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf */
template<class Type>
matrix<Type> logBeta(vector<Type> delta, matrix<Type> trMat, matrix<Type> lnProbs, int nbSteps) {
  int nbStates = trMat.cols();
  Type logbeta;
  matrix<Type> elnbeta(nbSteps,nbStates);
  matrix<Type> ltrMat(nbStates,nbStates);
  for(int j=0; j < nbStates; j++){
    elnbeta(nbSteps-1,j) = Type(0.0);
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
    }
  }
  for(int t=(nbSteps-2); t >= 0; t--){
    for(int i=0; i < nbStates; i++){
      logbeta = LOGZERO;
      for(int j=0; j < nbStates; j++){
        logbeta = elnsum(logbeta,elnproduct(ltrMat(i,j),elnproduct(lnProbs(t+1,j),elnbeta(t+1,j))));
      }
      elnbeta(t,i) = logbeta;
    }
  }
  return elnbeta;
}

/* Numerically stable state probabilities based on http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf */
template<class Type>
matrix<Type> stProbs(matrix<Type> elnalpha, matrix<Type> elnbeta, int nbSteps) {
  int nbStates = elnalpha.cols();
  Type normalizer;
  matrix<Type> gamma(nbSteps,nbStates);

  for(int t=0; t < nbSteps; t++){
    normalizer = LOGZERO;
    for(int i=0; i < nbStates; i++){
      gamma(t,i) = elnproduct(elnalpha(t,i),elnbeta(t,i));
      normalizer = elnsum(normalizer,gamma(t,i));
    }
    for(int i=0; i < nbStates; i++){
      gamma(t,i) = eexp(elnproduct(gamma(t,i),-normalizer));
    }
  }
  return gamma;
}

template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // DATA
    DATA_MATRIX(Y);	            //  (x, y) observations
    DATA_VECTOR(dt);         	//  time diff in some appropriate unit. this should contain dt for both interp and obs positions.
    DATA_VECTOR(mu0);        //  initial state mean
    DATA_SCALAR(V0);          // Prior variance of first location
    DATA_IVECTOR(isd);          //  indexes observations vs. interpolation points
    DATA_IVECTOR(obs_mod);      //  indicates which obs error model to be used
    DATA_INTEGER(nbStates);  // number of states
    DATA_INTEGER(nbSteps);  // number of predicted intervals 
    DATA_IVECTOR(nbObs);  // number of observations and predicted locations per time step
    DATA_INTEGER(secondObs); // time step for the 2nd observation
    
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
    
    // initial state and t.p.m. parameters
    PARAMETER_VECTOR(l_delta); // initial distribution (mlogit scale) - should be of length nbStates - 1
    PARAMETER_MATRIX(l_gamma); // t.p.m. (mlogit scale) - should be 1 x (nbStates-1) * nbStates
    
    // Transform parameters
    vector<Type> beta(nbStates);
    vector<Type> sigma2(nbStates);
    vector<Type> sigma(nbStates);
    for(int i=0; i < nbStates; i++){
      beta(i) = exp(log_beta(i));
      sigma2(i) = exp(2*log_sigma(i));
      sigma(i) = sqrt(sigma2(i));
    }
    Type psi = exp(l_psi);
    vector<Type> tau = exp(l_tau);
    Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
  
    vector<Type> delta(nbStates);
    matrix<Type> trMat(nbStates,nbStates);
    delta(0) = Type(1.0);
    for(int i=1; i < nbStates; i++){
      delta(i) = exp(l_delta(i-1));
    }
    delta /= delta.sum();
    int cpt = 0;
    for(int i=0;i<nbStates;i++){
      for(int j=0;j<nbStates;j++){
        if(i==j) {
          trMat(i,j) = Type(1.0);
          cpt++;
        } else trMat(i,j) = exp(l_gamma(0,i*nbStates+j-cpt));
      }
      trMat.row(i) /= trMat.row(i).sum();
    }
  
    int timeSteps = dt.size();
    
    /* Define likelihood */
    //parallel_accumulator<Type> jnll(this);
    matrix<Type> mu_bar_i(2,nbStates);
    vector<Type> sig_i(nbStates);
    Type nll_obs=0.0;
    
    /* Define likelihood */
    parallel_accumulator<Type> nll(this);
    //Type nll=0.0;
    
    matrix<Type> allProbs(nbSteps,nbStates);
    allProbs.setZero();
    vector<Type> delta0(nbStates);
    
    // CRW
  	// crwHMM PROCESS MODEL
    	
    for(int state = 0; state < nbStates; state++){

      //initial locs
      // time 0
      //delta0(state) = eexp(elnproduct(eln(delta(state)),dnorm(mu(0,0), mu0(0), V0,  1) + dnorm(mu(1,0), mu0(1), V0,  1)));
      delta0(state) = delta(state);
      
      //time 0 
      allProbs(0,state) += dnorm(mu(0,0), mu0(0), V0,  1) + dnorm(mu(1,0), mu0(1), V0,  1);
        
      //time 1
      sig_i(state) = dt(1)*sigma(state)/(Type(2.0)*beta(state));
      allProbs(secondObs,state) +=  dnorm(mu(0,1), mu(0,0), sig_i(state),  1) + dnorm(mu(1,1), mu(1,0), sig_i(state),  1); 
      
      int count = 1;
  	  for(int t=0; t < nbSteps; t++){
        for(int i = 0; i < nbObs[t]; i++) {
          count++;
          mu_bar_i.col(state) = (Type(1.0)+dt(count)/dt(count-1)-beta(state)*dt(count))*mu.col(count-1) + 
              (beta(state)*dt(count)-dt(count)/dt(count-1))*mu.col(count-2);
          sig_i(state) = dt(count)*sqrt(dt(count-1))*sigma(state);
          allProbs(t,state) += dnorm(mu(0,count), mu_bar_i(0,state), sig_i(state), 1) + 
            dnorm(mu(1,count), mu_bar_i(1,state), sig_i(state), 1);
        }
      }
    }

    REPORT(allProbs);

    // OBSERVATION MODEL
    // 2 x 2 covariance matrix for observations
    matrix<Type> cov_obs(2, 2);
    cov_obs.setZero();
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
      // forward algorithm
    Type hmmll = forward_alg(delta0, trMat, allProbs, nbSteps);
    REPORT(hmmll);
    nll += hmmll + nll_obs;

    ADREPORT(beta);
    ADREPORT(sigma);
    ADREPORT(rho_o);
    ADREPORT(tau);
    ADREPORT(psi);
  
    if(nbStates>1){
      vector<int> states = viterbi(delta0, trMat, allProbs, nbSteps);
      REPORT(states);
      matrix<Type> lalpha = logAlpha(delta0, trMat, allProbs, nbSteps);
      REPORT(lalpha);
      matrix<Type> lbeta = logBeta(delta0, trMat, allProbs, nbSteps);
      REPORT(lbeta);
      matrix<Type> stateProbs = stProbs(lalpha, lbeta, nbSteps); 
      REPORT(stateProbs);
      ADREPORT(delta);
      ADREPORT(trMat);
    }
    REPORT(nll_obs);
    
    return nll;

}

