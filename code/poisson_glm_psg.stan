// A poisson GLM
// I have also included additional components for the negative binomial
// variant which has one additional parameter
// these are currently commented out, but can be added back in
data {
  // Data declarations
  int N_obs;            // number of observations
  real Year[N_obs];     // vector of reals predictor values
  int Murres[N_obs];    // vector of integer counts
}

parameters {
  // Model parameter declarations
  real beta0;       // Defining intercept parameter as a real
  real beta1;       // Defining slope parameter as a real
  // real<lower=0> theta;    // Neg-bin version: theta = NB overdispersion parameter 
}

transformed parameters {
  // Parameter transformation
}

model {
  // prior definitions
  beta0 ~ normal(8, 4);       // Prior for intercept
  beta1 ~ normal(0, 0.8);     // Prior for slope
  // theta ~ normal(0,100);   // Neg-bin version: prior for theta
  
  // deterministic (i.e. the maths part)
  real lambda[N_obs];         // vector of means
    
  for(i in 1:N_obs){     // loop over observations
    // For each Count response, calculate lambda
    lambda[i] = exp(beta0 + (beta1 * Year[i]));   
  }
 
  // likelihood
  // Murres ~ neg_binomial_2(lambda, theta);     // Neg-bin version 
  Murres ~ poisson(lambda);        // comment this line out if you want the NB version
}

generated quantities{
  real c1998;
  c1998 = exp(beta0 + (beta1 * (1998-2005)));
}
