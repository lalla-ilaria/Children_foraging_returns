data{
	int N;      //number of children
	int M;      //number of trip
	int Q;      //number of questions
	int ID_i[M];//index of forager/return
	int S[M];   //success of trap
	real L[M];  //length of trip
	int Y[N,Q]; //answers
	real H[N];  //individual height
	real G[N];  //individual grip
	real A[N];  //individual age
	}
parameters{
  vector<lower=0>[Q] a;
  vector<lower=0>[Q] b;
  vector<lower=0>[N] K;
  vector [N] iota; //individual level random effect
  real<lower=0> sigma_i;
  real<lower=0> alpha;
  real<lower=0> beta; //age effect
  real<lower=0> gamma; //age elasticity
  real<lower=0> zeta_k; //knowledge elasticity
  real<lower=0> eta_h; //height elasticity
  real<lower=0> theta_g; //grip elasticity
  real xi; //exponent for length trip
}
transformed parameters{
  vector [N] phi; //individual total effect
  vector [M] psi; //trip effect
  for(i in 1:N) phi[i]  = exp (iota[i] * sigma_i) * ( 
                         (1-exp(-beta * A[i]  )) ^ gamma * 
                          1 + K[i] ^ zeta_k * 
                          H[i] ^ eta_h * 
                          G[i] ^ theta_g);
  for(i in 1:M) psi[i] =  L[i] ^ xi;

}
model{
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma~ lognormal(0, 1);
  zeta_k~ lognormal(0, 1);
  eta_h ~ lognormal(0, 1);
  theta_g ~ lognormal(0, 1);
  xi ~ normal(0, 1);
 	a ~ lognormal(0, 0.5); //value constrained above zero
	for(i in 1:Q) b[i] ~ normal(3.3,0.5);
	for(i in 1:N) K[i] ~ normal(3.3,0.5);
	for ( i in 1:N ) {
		for (j in 1:Q ) {
			real p_irt = inv_logit(a[j]*(K[i]-b[j]));
			Y[i,j] ~ bernoulli( p_irt );
		}
	}
  for ( i in 1:M ) {
         real p = 1 - exp (- alpha * phi[ID_i[i]] * psi[i]);
         S[i] ~ bernoulli(p); 
      }
}
