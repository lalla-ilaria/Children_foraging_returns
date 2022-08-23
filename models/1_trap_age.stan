data{
	int N;      //number of children
	int M;      //number of trip
	array[M] int ID_i;//id of forager/return
	array[M] int success;  //returns
	array[M] real duration;  //length of trip
	array[N] real age;
	}
parameters{
  vector [N] iota;
  real<lower=0> sigma_i;
  real<lower=0> alpha;
  real<lower=0> beta; //age effect
  real<lower=0> gamma; //age elasticity
	real xi; //exponent for length trip
}
transformed parameters{
  vector [N] phi;
  vector [M] psi;
  for(i in 1:N) phi[i]  = iota[i] * sigma_i + 
                          gamma * log(1-exp(-beta * age[i]  ));
  for(i in 1:M) psi[i] =  xi * log(duration[i]);

}
model{
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma ~lognormal(0, 1);
  xi ~ normal(0, 1);
  for ( i in 1:M ) {
         real p = 1 - exp (- alpha * exp(phi[ID_i[i]]) * exp(psi[i]));
         success[i] ~ bernoulli(p); 
      }
}
