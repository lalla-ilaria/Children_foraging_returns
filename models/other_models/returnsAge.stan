data{
	int M;      //number of trip
	real R[M];  //returns
	real L[M];  //length of trip
	real K[M];  //individual knowledge
	real B[M];  //individual body
	real A[M];
	}
parameters{
  real alpha_k;
  real k_a;
  real sigma_k;
  real alpha_b;
  real b_a;
  real sigma_b;
  real alpha_r;
  real<lower=0> r_a; //age effect
  real<lower=0> r_l; //exponent for length of trip - amount per trip
	real<lower=0> r_k; //exponent for knowledge - amount per trip
	real<lower=0> r_b; //exponent for skill - amount per trip
	real<lower=0> sigma_r;
}
model{
  alpha_k ~ normal(0,1);
  k_a ~ normal(0,1);
  sigma_k ~ exponential(1);
  alpha_b ~ normal(0,1);
  b_a ~ normal(0,1);
  sigma_b ~ exponential(1);
  alpha_r ~ normal(0,1);
  r_a ~ lognormal(0,2);
  r_l ~ lognormal(0,2);
  r_k ~ lognormal(0,2);
  r_b ~ lognormal(0,2);
  sigma_r ~ exponential(1);
  for ( i in 1:M ) {
    real mu_k = alpha_k +  A[i] ^ k_a ;
    real mu_b = alpha_b +  A[i] ^ b_a ;
    real mu_r = alpha_r + 
             r_a * log(A[i]) +
             r_l * log(L[i]) + r_k * log(K[i]) + r_b * log(B[i]) ;
    K[i] ~ normal ( mu_k , sigma_k) ;
    B[i] ~ normal ( mu_b , sigma_b) ;
		R[i] ~ lognormal( exp(mu_r) , sigma_r ); 
      }
}
