data{
	int N;
	int M;      //number of trip
	real R[M];  //returns
	real A[M];  //individual age
	int ID[M];
	}
parameters{
  real alpha;
  vector[N] alpha_i;
  real<lower=0> r_a; //age effect
	real<lower=0> sigma;
	real<lower=0> sigma_a;
}
model{
  vector [M] m;
  alpha ~ normal (0, 0.1); 
  r_a ~ normal(0,1);
  sigma ~ exponential(1);
  sigma_a ~ exponential(5);
  alpha_i ~ lognormal (exp(alpha), sigma_a); 
  for ( i in 1:M ) {
    m[i] = alpha_i[ID[i]] ;
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