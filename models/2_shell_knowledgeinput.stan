data{
	//foraging data
	int N;            //number of individual in foraging data
	int M;            //number of trip
	int ID_i[M];      //index of forager/return
	int ID_k[N];      //index of forager/return for knowledge data
	real returns[M];  //returns
	real duration[M]; //length of trip
	real tide[M];     //height of tide
	real age[N];       //individual age
	int age_int[N];       //individual age
	int sex[N];       //individual sex
	// real height[N];     //individual height
	// real grip[N];     //individual grip
	int knowledge_impute[N]; //whether knowledge data are available or not
	
	//knowledge data
	int W;            //number of individuals in knowledge data
	int Q;            //number of questions
	int O;            //n ages prior drichlet
	int age_irt[W];   //age of individuals
	int sex_irt[W];   //sex of indiviudals
  int answers[W,Q]; //answers to freelist
	vector[O-1] prior_dirichlet; //prior drichlet
	}
parameters{
  //foraging parameters
  vector [N] iota;       //individual level random effect
  real<lower=0> sigma_i; //sd of iota
  real<lower=0> alpha;   //scaling intercept
  real<lower=0> beta;    //age effect
  real<lower=0> gamma;   //age elasticity
  real<lower=0> zeta_k;  //knowledge elasticity
  // real<lower=0> eta_h; //height elasticity
  // real<lower=0> theta_g; //grip elasticity
  real xi;               //exponent for length trip
  real tau;              //coefficient of tide
	real<lower=0> sigma;   //sd of lognormal
	
	//knowledge parameters
	real omega;
	vector[W] iota_irt;       // individual intercepts on knowledge
  real<lower=0> sigma_irt;  //sd for iota_irt
  vector<lower=0>[2] ro_age; // coefficient relating age to knowledge, one per sex
  vector<lower=0>[Q] a;     //discrimination of questions
  vector<lower=0>[Q] b;     //difficulty of questions
  simplex[O-1] delta ;       //age specific effects
}
transformed parameters{
  //additional parameters
  vector<lower=0>[W] knowledge;
  vector<lower=0>[N] knowledge_merge;
  vector[O] delta_j;
  vector[N] phi;                   //individual total effect
  vector[M] psi;                   //trip effect
  delta_j  = append_row(0, delta);
  
  //knowledge estimation
  for ( i in 1:W )   knowledge[i] = omega + iota_irt[i] * sigma_irt +                    //individual interecepts -absorbs residual variation   
                                    ro_age[sex_irt[i]] * sum (delta_j[ 1 : age_irt[i]] );//effect of age - sex specific
  //impute knowledge from irt if missing data
  for(i in 1:N) {
    if (knowledge_impute[i] == 0) knowledge_merge[i] = knowledge[ID_k[i]]; 
    if (knowledge_impute[i] == 1) knowledge_merge[i] = omega + ro_age[sex[i]] * sum (delta_j[ 1 : age_int[i]] ); 
  }//i
  
  //foraging estimation
  for(i in 1:N) phi[i]  = exp (iota[i] * sigma_i) * ( 
                         (1-exp(-beta * age[i]  )) ^ gamma * 
                          knowledge_merge[i] ^ zeta_k );//* 
                          // height[i] ^ eta_h * 
                          // grip[i] ^ theta_g);
  for(i in 1:M) psi[i] =  duration[i] ^ xi * 
                          exp( tide[i] * tau);//height of tide
}
model{
  //priors for foraging
  for(i in 1:N) knowledge_merge[i] ~ normal(4,0.5) ;
 	iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma~ lognormal(0, 1);
  zeta_k~ lognormal(0, 1);
  // eta_h ~ lognormal(0, 1);
  // theta_g ~ lognormal(0, 1);
  xi ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ exponential(1);
  //priors for knowledge
  omega ~ normal (0,2);
  for(i in 1:W) knowledge[i] ~ normal(4,0.5) ;
 	a ~ lognormal(0, 0.5); //value constrained above zero
	for(i in 1:Q) b[i] ~ normal(4,0.5);
	iota_irt ~ normal(0,0.1);
  ro_age ~ lognormal( 0 , 1 );
  delta ~ dirichlet( prior_dirichlet );
  sigma_irt ~ exponential(1);

	//irt across all people to estimage age and sex specific knowledge
	for ( i in 1:W ) {
		for (j in 1:Q ) {
			real p_irt = inv_logit(a[j]*(knowledge[i]-b[j]));
			answers[i,j] ~ bernoulli( p_irt );
		}
	}
  for ( i in 1:M ) {
         real m = log( alpha * phi[ID_i[i]] * psi[i]);
         returns[i] ~ lognormal( m , sigma ); 
      }
}
