data{
	int N;      //number of children
	int M;      //number of trip
	array[M] int ID_i;//id of forager/return
	array[M] real returns;  //returns
	array[M] real duration;  //length of trip
	array[N] real age;
	array[M] real tide;//height of tide
	}
parameters{
  vector [N] iota;
  real<lower=0> sigma_i;
  real<lower=0> alpha;
  real<lower=0> beta; //age effect
  real<lower=0> gamma; //age elasticity
  real xi; //exponent for length trip
  real tau;//coefficient of tide
	real<lower=0> sigma;
}
transformed parameters{
  vector [N] phi;
  vector [M] psi;
  for(i in 1:N) phi[i]  = iota[i] * sigma_i +  
                          gamma * log(1-exp(-beta * age[i]  )) ;
  for(i in 1:M) psi[i] =  xi * log(duration[i]) + 
                          tau * tide[i] ;//height of tide
;

}
model{
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma~ lognormal(1, 1);
  xi ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ exponential(1);
  for ( i in 1:M ) {
         real m = log( alpha) + phi[ID_i[i]] + psi[i];
         returns[i] ~ lognormal( m , sigma ); 
      }
}
