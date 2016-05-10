// Based on Supplementary Web Appendix 1 to:
// Malcolm J. Price, A. E. Ades, Daniela De Angelis, Nicky J. Welton, John Macleod, Kate Soldan, Katy Turner, Ian Simms and Paddy J. Horner
// Mixture-of-exponentials models to explain heterogeneity in studies of the duration of Chlamydia trachomatis infection.
// Statistics in Medicine 32:1547–1560 (2013)

data {
  int<lower=0> studnum; // number of studies 
  int<lower=0> studnum_bytype[3]; // number of studies 
  int<lower=0> studobs[studnum]; // number of observations (time periods) in each study
  int<lower=0> cumobs[studnum]; // cumulative number of observations at the start of each study
  int<lower=0> Nobs; // total number of observations (=sum(studobs))
  int<lower=0> r[Nobs]; // number who cleared infection at each time point 
  int<lower=0> n[Nobs]; // number tested at each time point 
  real<lower=0> t[Nobs]; // estimated mean follow-up 
  int<lower=0> seind[Nobs]; // did the study use culture as opposed to NAATs?
  real<lower=0> T[Nobs]; // already followed for...
}

parameters {
  real<lower=0,upper=1> p1; // probability of clearing fast
  real<lower=0> lambda_slow; // fast and slow clearance rates
  real<lower=0,upper=1> psi; // sensitivity of culture given initial positive culture
}

transformed parameters {
	
  real<lower=0,upper=1> w1; //  expected proportion of fast clearers in a prevalent population
  real theta[Nobs]; // weighted average of clearance probabilities for two classes
  real<lower=0> lambda[2];
  real<lower=0,upper=1> pk; // temporary variable, containing probability of this k
  
  // rates
  lambda[1] <- 120;
  lambda[2] <- lambda_slow;
  
  // proportion of participants in each study time point expected to have recovered.
  // clinic studies
  for (i in 1:studnum_bytype[1]) { // should find a way of using input data to state study categories
    for (j in 1:studobs[i]) {
      theta[cumobs[i]+j] <-  0;
      for(k in 0:n[cumobs[i]+j]){
        theta[cumobs[i]+j] <- theta[cumobs[i]+j] + exp( binomial_coefficient_log(n[cumobs[i]+j],k) ) * p1^k * (1-p1)^(n[cumobs[i]+j]-k) *  ((k / (1.0*n[cumobs[i]+j])) * (1 - exp(-lambda[1] * t[cumobs[i]+j])) + (1 - (k / (1.0*n[cumobs[i]+j]))) * (1 - exp(-lambda[2] * t[cumobs[i]+j])));
        }
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
   
  // left-truncated; single observation
  w1 <- (p1 / lambda[1]) / (p1 / lambda[1] + (1 - p1) / lambda[2]);
  for (i in (studnum_bytype[1]+1):(studnum_bytype[1]+studnum_bytype[2])) { // should find a way of using input data to state study categories
    for (j in 1:studobs[i]) {
      theta[cumobs[i]+j] <- 0;
      for(k in 0:n[cumobs[i]+j]){
        theta[cumobs[i]+j] <- theta[cumobs[i]+j] + exp( binomial_coefficient_log(n[cumobs[i]+j],k) ) * w1^k * (1-w1)^(n[cumobs[i]+j]-k) *  ((k / (1.0*n[cumobs[i]+j])) * (1 - exp(-lambda[1] * t[cumobs[i]+j])) + (1 - (k / (1.0*n[cumobs[i]+j]))) * (1 - exp(-lambda[2] * t[cumobs[i]+j])));
        }
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);
      }
    }
    
  // left-truncated; repeat observations
  for (i in (studnum_bytype[1]+studnum_bytype[2]+1):studnum) {
    for (j in 1:studobs[i]) {
        theta[cumobs[i]+j] <- 0;
        for (k in 0:n[cumobs[i]+j]){
          
          pk <- exp( binomial_coefficient_log(n[cumobs[i]+j],k) ) * w1^k * (1-w1)^(n[cumobs[i]+j]-k); // temporary variable, containing probability of this k

          theta[cumobs[i]+j] <- theta[cumobs[i]+j] + pk * (k / (1.0*n[cumobs[i]+1])) * exp(-lambda[1]*T[cumobs[i]+j]) * (1 - exp(-lambda[1]*t[cumobs[i]+j])) / ( (k / (1.0*n[cumobs[i]+1])) * exp(-lambda[1]*T[cumobs[i]+j]) + (1 - k / (1.0*n[cumobs[i]+1])) * exp(-lambda[2]*T[cumobs[i]+j]) );   // fast category
          theta[cumobs[i]+j] <- theta[cumobs[i]+j] + pk * (1 - k / (1.0*n[cumobs[i]+1])) * exp(-lambda[2]*T[cumobs[i]+j]) * (1 - exp(-lambda[2]*t[cumobs[i]+j])) / ( (k / (1.0*n[cumobs[i]+1])) * exp(-lambda[1]*T[cumobs[i]+j]) + (1 - k / (1.0*n[cumobs[i]+1])) * exp(-lambda[2]*T[cumobs[i]+j]) ) ;  // slow category
          
          }
      if(seind[cumobs[i]+j])
        theta[cumobs[i]+j] <- 1 + psi*(theta[cumobs[i]+j] - 1);    
      }
    }
  
    
  }

model {

  // priors
  p1 ~ beta(1,1); // proportion in fast-clearing component 
  lambda_slow ~ exponential(0.001); // slow clearance rate
  psi ~ beta(78,8); // sensitivity of culture given initial positive culture
	
  // Class proportions
  
//  // t=0 studies
//  for (i in 1:4) {	
//    for (j in 1:studobs[i]) {	
//      z_int[i,j] <- trunc(z[i,j]);
//      z_int[i,j] ~ binomial(n[i,j], p1); // start at t=0
//      }
//    }

//  // Left-truncated studies		
//  for (i in 5:studnum) {	
//    for (j in 1:studobs[i]) {	
//      z_int[i,j] <- floor(z[i,j]);
//      z_int[i,j] ~ binomial(n[i,j], w1);  
//      }
//    }
	

  // likelihood
  for (i in 1:studnum) {
    for (j in 1:studobs[i]) {
      print(i);
      print(j);
      r[cumobs[i]+j] ~ binomial(n[cumobs[i]+j], theta[cumobs[i]+j]);
      }
    }

  }

generated quantities {

//  // deviance
//  for (i in 1:studnum) {
//    for (j in 1:studobs[i]) {
//    
//      dev[i,j] <- 2 * (r[i,j] * log(r[i,j] / (theta[i,j] * n[i,j])) + (n[i,j] - r[i,j]) * log((n[i,j] - r[i,j]) / (n[i,j] - (n[i,j] * theta[i,j]))))	
//      
//      }
//    dev.stud[i] <- sum(dev[i,1:studobs[i]])
//    }
//  sumdev <- sum(dev.stud[])
//	
//  // summary statistics
//  dur <- 1 / lambda[2]

//  // Predicted values for Forest plot
//  for (i in 1:studnum) {
//    for (j in 1:studobs[i]) {	
//      stud.lambdaexpect[i,j] <- -log(1 - theta[i,j]) / t[i,j]
//      stud.dur.expect[i,j] <- 1 / stud.lambdaexpect[i,j]	
//      }
//    }
  
}

