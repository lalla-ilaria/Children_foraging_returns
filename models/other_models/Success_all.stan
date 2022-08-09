data{
	int N;      //number of children
	int M;      //number of trip
	int ID_i[M];//id of forager/return
	int S[M];  //returns
	real L[M];  //length of trip
	real K[N];  //individual knowledge
	real B[N];  //individual knowledge
	real A[N];
	}
parameters{
  vector [N] iota;
  real<lower=0> sigma_i;
  real<lower=0> alpha;
  real<lower=0> beta; //age effect
  real<lower=0> gamma; //age elasticity
  real zeta_k; //knowledge elasticity
  real eta_b; //knowledge elasticity
	real xi; //exponent for length trip
}
transformed parameters{
  vector [N] phi;
  vector [M] psi;
  for(i in 1:N) phi[i]  = exp(iota[i] * sigma_i) * 
                          (1-exp(-beta * A[i]  )) ^ gamma * 
                          K[i] ^ zeta_k *
                          B[i] ^ eta_b;
  for(i in 1:M) psi[i] =  L[i] ^ xi;

}
model{
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma ~lognormal(0, 1);
  zeta_k ~ normal(0, 1);
  eta_b ~ normal(0, 1);
  xi ~ normal(0, 1);
  for ( i in 1:M ) {
         real p = 1 - exp (- alpha * phi[ID_i[i]] * psi[i]);
         S[i] ~ bernoulli(p); 
      }
}

  // phi[i,] <- log(1-exp(-beta_a[i] * AGE  )) * gamma_a[i]
  // psi[i] <-   lambda[i] * log (L[3])
  // R <- exp ( alpha[i] + phi[i,] + psi[i] + ((sigma[i]^2) /2))

