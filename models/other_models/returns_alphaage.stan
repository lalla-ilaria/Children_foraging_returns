data{
	int N;
	int M;      //number of trip
	real R[M];  //returns
	real L[M];  //length of trip
	real K[M];  //individual knowledge
	real B[M];  //individual body
	real A[M];  //individual age
	int ID[M];
	}
parameters{
  real alpha;
  vector[N] alpha_i;
  real<lower=0> r_a; //age effect
  real<lower=0> r_l; //exponent for length of trip - amount per trip
	real<lower=0> r_k; //exponent for knowledge - amount per trip
	real<lower=0> r_b; //exponent for skill - amount per trip
	real<lower=0> sigma;
	real<lower=0> sigma_a;
}
model{
//  alpha ~ normal (0,1);
  vector [M] m;
  vector [M] mu_a;
  alpha ~ normal (0, 1); 
  r_a ~ normal(0,1);
  r_l ~ lognormal(0,2);
  r_k ~ lognormal(0,2);
  r_b ~ lognormal(0,2);
  sigma ~ exponential(1);
  sigma_a ~ exponential(1);
  for ( i in 1:M ) {
    mu_a[i] = alpha + r_a * A[i];
    alpha_i[ID[i]] ~ normal (mu_a[i], sigma_a); 
    m[i] = alpha_i[ID[i]] + 
             r_l * log(L[i]) + r_k * log(K[i]) + r_b * log(B[i]) ;
		R[i] ~ lognormal( exp(m) , sigma ); 
      }
}
// 
//     vector[48] p;
//     sigma ~ exponential( 1 );
//     a_bar ~ normal( 0 , 1.5 );
//     a ~ normal( a_bar , sigma );
//     for ( i in 1:48 ) {
//         p[i] = a[tank[i]];
//         p[i] = inv_logit(p[i]);
//     }
//     S ~ binomial( N , p );