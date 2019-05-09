#define TMB_LIB_INIT R_init_mypkg
#include <TMB.hpp>
#include <cmath>

/*
	Random walk with velocities as a random walk.
	Kind of a correlated RW model - with code for KF and LS argos errors

	Authors:
		Ian Jonsen 	(ian.jonsen@mq.edu.au)
		Toby Patterson (Toby.Patterson@csiro.au)
*/

using namespace density;
using std::sqrt;

/* implement the vector - matrix product */
template<class Type>
vector<Type> multvecmat(array<Type>  A, matrix<Type>  B) {
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
  else C = sqrt(-2.0);
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
    C = sqrt(-2.0);
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
  Type nan = sqrt(-2.0);
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  vector<Type> lalpha(nbStates);
  vector<Type> lnewalpha(nbStates);
  Type sumalpha = nan;
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = nan;
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
    sumalpha = nan;
    for(int j=0; j < nbStates; j++){
      logalpha = nan;
      for(int i=0; i < nbStates; i++){
        logalpha = elnsum(logalpha,elnproduct(lalpha(i),ltrMat(i,j)));
      }
      lnewalpha(j) = elnproduct(logalpha,lnProbs(t,j));
      sumalpha  = elnsum(sumalpha,lnewalpha(j));
    }
    jnll -= sumalpha;
    lnewalpha -= sumalpha;
    lalpha = lnewalpha;
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
  Type nan = sqrt(-2.0);
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = nan;
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
  Type nan = sqrt(-2.0);
  matrix<Type> ltrMat(nbStates,nbStates);
  vector<Type> ldeltaG(nbStates);
  for(int j=0; j < nbStates; j++){
    ldeltaG(j) = nan;
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
      ldeltaG(j) = elnsum(ldeltaG(j),elnproduct(eln(delta(i)),ltrMat(i,j)));
    }
    elnalpha(0,j) = elnproduct(ldeltaG(j),lnProbs(0,j));
  }
  for(int t=1; t < nbSteps; t++){
    for(int j=0; j < nbStates; j++){
      logalpha = nan;
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
  Type nan = sqrt(-2.0);
  matrix<Type> ltrMat(nbStates,nbStates);
  for(int j=0; j < nbStates; j++){
    elnbeta(nbSteps-1,j) = Type(0.0);
    for(int i=0; i < nbStates; i++){
      ltrMat(i,j) = eln(trMat(i,j));
    }
  }
  for(int t=(nbSteps-2); t >= 0; t--){
    for(int i=0; i < nbStates; i++){
      logbeta = nan;
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
  Type nan = sqrt(-2.0);

  for(int t=0; t < nbSteps; t++){
    normalizer = nan;
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
    DATA_VECTOR(state0);        //  initial state
    DATA_IVECTOR(isd);          //  indexes observations vs. interpolation points
    DATA_IVECTOR(obs_mod);      //  indicates which obs error model to be used
    DATA_INTEGER(proc_mod);		//	indicates which process model to be used: RW or CRW
    DATA_INTEGER(nbStates);  // number of states
    DATA_INTEGER(nbSteps);  // number of predicted intervals 
    DATA_IVECTOR(nbObs);  // number of observations and predicted locations per time step

    // for KF observation model
    DATA_VECTOR(m);             //  m is the semi-minor axis length
    DATA_VECTOR(M);             //  M is the semi-major axis length
    DATA_VECTOR(c);             //  c is the orientation of the error ellipse
    // for LS observation model
    DATA_MATRIX(K);                 // error weighting factors for LS obs model

    // PROCESS PARAMETERS
    // for RW
    PARAMETER_MATRIX(l_sigma);    //  Innovation variance (link scale)
    PARAMETER_VECTOR(l_rho_p);           //  Innovation correlation (link scale)
    PARAMETER_MATRIX(X);          //  Predicted locations TP - length(X) should be same as length(dt) - i.e. both interp & obs pos.
    // for CRW
    PARAMETER_VECTOR(logD);				// 1-d Diffusion coefficient
    // random variables
    PARAMETER_MATRIX(mu); /* State location */
    PARAMETER_MATRIX(v); /* state velocities */

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
	  matrix<Type> sigma(nbStates,2);
	  vector<Type> delta(nbStates);
	  matrix<Type> trMat(nbStates,nbStates);
	  if(nbStates>1){
  	  delta(0) = Type(1.0);
	    for(int i=1; i < nbStates; i++){
	      delta(i) = exp(l_delta(i-1));
	    }
  	  delta /= (delta.sum());
  	  vector<Type> rowSums(nbStates);
  	  rowSums.setZero();
  	  int cpt = 0;
  	  for(int i=0;i<nbStates;i++){
        for(int j=0;j<nbStates;j++){
          if(i==j) {
            trMat(i,j) = Type(1.0);
            cpt++;
          } else trMat(i,j) = exp(l_gamma(0,i*nbStates+j-cpt));
          // keep track of row sums, to normalize in the end
          rowSums(i) += trMat(i,j);
        }
        for(int j=0;j<nbStates;j++){
          trMat(i,j) /= rowSums(i);
        }
  	  }
	  } else {
	    delta(0) = Type(1.0);
	    trMat(0,0) = Type(1.0);
	  }
	  
	  for(int state = 0; state < nbStates; state++){
	    sigma(state,0) =  exp(l_sigma(state,0));
	    sigma(state,1) =  exp(l_sigma(state,1));
	  }
    vector<Type> rho_p = Type(2.0) / (Type(1.0) + exp(-l_rho_p)) - Type(1.0);
    vector<Type> tau = exp(l_tau);
    Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);
    Type psi = exp(l_psi);
    vector<Type> D = exp(logD);

    int timeSteps = dt.size();

	/* Define likelihood */
    parallel_accumulator<Type> jnll(this);
    //Type jnll = 0.0;
    matrix<Type> allProbs(nbSteps,nbStates);
    allProbs.setZero();
    matrix<Type> probs(timeSteps-1,nbStates);
    Type tiny = 1e-5;

    if(proc_mod == 0) {
    	// RW
    	// 2 x 2 covariance matrix for innovations
    	matrix<Type> cov(2, 2);
    	matrix<Type> cov_dt(2, 2);            // tmp matrix for dt * cov calcs withn process loop
		  MVNORM_t<Type> nll_proc(cov);

    	// RW PROCESS MODEL
      	
      for(int state = 0; state < nbStates; state++){
        
        cov(0, 0) = sigma(state,0) * sigma(state,0);
        cov(0, 1) = rho_p(state) * sigma(state,0) * sigma(state,1);
        cov(1, 0) = cov(0, 1);
        cov(1, 1) = sigma(state,1) * sigma(state,1);
        
      	int count = 1;
    	  for(int t=0; t < nbSteps; t++){
          for(int i = 0; i < nbObs[t]; i++) {
        		cov_dt = dt(count) * dt(count) * cov;
        		nll_proc.setSigma(cov_dt);
        		probs(count-1,state) = -nll_proc(X.row(count) - X.row(count - 1)); // Process likelihood
        		allProbs(t,state) += probs(count-1,state);
            count += 1;
        	}
      	}
      }
    	//REPORT(cov);
    	//REPORT(X);
    } else if(proc_mod == 1){
    	// CRW
		// Setup object for evaluating multivariate normal likelihood
    	matrix<Type> cov(4,4);
      MVNORM_t<Type> nll_proc(cov);

    	// CRW PROCESS MODEL
    	vector<Type> x_t(4);
    	
    	// loop over 2 coords and update nll of start location and velocities.
    	for(int i = 0; i < 2; i++) {
    	  jnll -= dnorm(mu(i,0), state0(i), tiny, true);
    	  jnll -= dnorm(v(i,0), state0(i+2), tiny, true);
    	}
    	
	    //cov.setZero();
	    //cov(0,0) = tiny;
	    //cov(1,1) = tiny;
	    //cov(2,2) = 2 * D(state) * dt(0);
    	//cov(3,3) = 2 * D(state) * dt(0);

    	// loop over 2 coords and update nll of start location and velocities.
    	//for(int i = 0; i < 2; i++) {
    	//  jnll -= dnorm(mu(i,0), state0(i), tiny, true);
    	//  jnll -= dnorm(v(i,0), state0(i+2), tiny, true);
    	//}
    	
    	for(int state=0; state < nbStates; state++){
    	  
      	int count = 1;
    	  for(int t=0; t < nbSteps; t++){
    	    
        	for(int i = 0; i < nbObs[t]; i++) {
        	  
        	  // location innovations
        	  x_t(0) = mu(0,count) - (mu(0,count-1) + (v(0,count) * dt(count)));
        	  x_t(1) = mu(1,count) - (mu(1,count-1) + (v(1,count) * dt(count)));
        	  
        	  // velocity innovations
        	  x_t(2) = (v(0,count) - v(0,count-1)); // /dt(count);
        	  x_t(3) = (v(1,count) - v(1,count-1)); // /dt(count);
        	  
        		// process cov at time t
        		cov.setZero();
        		cov(0,0) = tiny;
        		cov(1,1) = tiny;
        		cov(2,2) = 2 * D(state) * dt(count);
        		cov(3,3) = 2 * D(state) * dt(count);
        		nll_proc.setSigma(cov);
        		probs(count-1,state) =  -nll_proc(x_t); // Process likelihood
        		allProbs(t,state) += probs(count-1,state);
        	  count += 1;
          }
    	  }
    	}
    	//REPORT(x_t);
    }
    
    REPORT(probs);
    REPORT(allProbs);

    // OBSERVATION MODEL
    // 2 x 2 covariance matrix for observations
    matrix<Type> cov_obs(2, 2);
    MVNORM_t<Type> nll_obs; // Multivariate Normal for observations

    for(int i = 0; i < timeSteps; ++i) {
      if(isd(i) == 1) {
        if(obs_mod(i) == 0) {
          // Argos Least Squares observations
          Type s = tau(0) * K(i,0);
          Type q = tau(1) * K(i,1);
          cov_obs(0,0) = s * s;
          cov_obs(1,1) = q * q;
          cov_obs(0,1) = s * q * rho_o;
          cov_obs(1,0) = cov_obs(0,1);
        } else if(obs_mod(i) == 1) {
          // Argos Kalman Filter (or Kalman Smoothed) observations
          double z = sqrt(2.);
          Type s2c = sin(c(i)) * sin(c(i));
          Type c2c = cos(c(i)) * cos(c(i));
          Type M2  = (M(i) / z) * (M(i) / z);
          Type m2 = (m(i) * psi / z) * (m(i) * psi / z);
          cov_obs(0,0) = (M2 * s2c + m2 * c2c);
          cov_obs(1,1) = (M2 * c2c + m2 * s2c);
          cov_obs(0,1) = (0.5 * (M(i) * M(i) - (m(i) * psi * m(i) * psi))) * cos(c(i)) * sin(c(i));
          cov_obs(1,0) = cov_obs(0,1);
        }
        nll_obs.setSigma(cov_obs);   // set up i-th obs cov matrix
        if(proc_mod == 0) {
          jnll += nll_obs(Y.row(i) - X.row(i));   // RW innovations
        } else if(proc_mod == 1) {
          jnll += nll_obs(Y.col(i) - mu.col(i));   // CRW innovations
        }

      }
    }
    
    // forward algorithm
    Type hmmll = forward_alg(delta, trMat, allProbs, nbSteps);
    REPORT(hmmll);
    jnll += hmmll;

    if(proc_mod == 0) {
      ADREPORT(rho_p);
      ADREPORT(sigma);
    } else if(proc_mod == 1) {
      ADREPORT(D);
    }
    ADREPORT(rho_o);
    ADREPORT(tau);
    ADREPORT(psi);
    if(nbStates>1){
      vector<int> states = viterbi(delta, trMat, allProbs, nbSteps);
      REPORT(states);
      matrix<Type> lalpha = logAlpha(delta, trMat, allProbs, nbSteps);
      REPORT(lalpha);
      matrix<Type> lbeta = logBeta(delta, trMat, allProbs, nbSteps);
      REPORT(lbeta);
      matrix<Type> stateProbs = stProbs(lalpha, lbeta, nbSteps); 
      REPORT(stateProbs);
      ADREPORT(delta);
      ADREPORT(trMat);
    }
    return jnll;
}

