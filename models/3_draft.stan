data{
	//individual level data
	int N;            //total number of individuals
	real age[N];       //individual age
	int sex[N];       //individual sex
	
	//data availability
	int has_foraging[N];
	int has_knowledge[N];
	
	// foraging trip data
	int M;            //number of trip
	int ID_i[M];      //index of forager/return
	real returns[M];  //returns
	real duration[M]; //length of trip
	real tide[M];     //height of tide
	
	//knowledge data
	int Q;            //number of questions
	int answers[N,Q]; //answers to freelist
}//data
	
parameters{
  //foraging parameters
  vector [N] iota;       //individual level random effect
  real<lower=0> sigma_i; //sd of iota
  real<lower=0> alpha;   //scaling intercept
  real<lower=0> beta;    //age effect
  real<lower=0> gamma;   //age elasticity
  real<lower=0> zeta_k;  //knowledge elasticity
  real xi;               //exponent for length trip
  real tau;              //coefficient of tide
	real<lower=0> sigma;   //sd of lognormal
	
	//knowledge parameters
	vector[N] knowledge;
	real<lower=0> sigma_k;  //sd for knowledge estimation
  real omega;
	vector[N] iota_irt;       //individual intercepts on knowledge
  real<lower=0> sigma_irt;  //sd for iota_irt
  vector<lower=0>[2] chi;   //max effect of age, sex specific
  vector<lower=0>[2] ro_age;//sex specific age speed
  vector<lower=0>[Q] a;     //discrimination of questions
  vector<lower=0>[Q] b;     //difficulty of questions
}//parameters

transformed parameters{
  //additional parameters
  vector[N] knowledge_st;
  vector[N] phi;                   //individual total effect
  vector[M] psi;                   //trip effect
  
  //standardize knowledge to fit it to foraging model
  for (i in 1:N) {
  knowledge_st[i] = (knowledge[i]-min(knowledge))/ mean ((knowledge - min(knowledge)));
  }
  
  //estimate individual and trip level parameters
  for(i in 1:N) {
    phi[i]  = exp (iota[i] * sigma_i) * (
                         (1-exp(-beta * age[i]  )) ^ gamma *
                          knowledge_st[i] ^ zeta_k );
    }
  for(i in 1:M) {
    psi[i] =  duration[i] ^ xi *
                          exp( tide[i] * tau);//height of tide
    }
}//transformed parameters
model{
  //priors for foraging
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ lognormal(0, 1);
  gamma~ lognormal(0, 1);
  zeta_k~ lognormal(0, 1);
  xi ~ normal(0, 1);
  tau ~ normal(0, 1);
  sigma ~ exponential(1);
  
  //priors for knowledge
  sigma_k ~ exponential(1);
  omega ~ normal (0,2);
  a ~ exponential(1); 
	for(i in 1:Q) b[i] ~ normal(1,0.5);
	iota_irt ~ normal(0,1);
  chi ~ exponential(1);
  ro_age ~ exponential(1);
  sigma_irt ~ exponential(1);

	//irt across knowledge people to estimage age and sex specific knowledge
	for ( i in 1:N ) {
		for (j in 1:Q ) {
		  if (has_knowledge[i] == 1){
		    real p_irt = inv_logit(a[j]*(knowledge[i]-b[j]));
			  answers[i,j] ~ bernoulli( p_irt );
		  }
		}
	}
  //estimate knowledge for all individuals
  for ( i in 1:N ){
        real m = omega + 
                 iota_irt[i] +  
                 chi[sex[i]] * ( 1 - exp(-ro_age[sex[i]] * age[i]));
        knowledge[i] ~ normal(m, sigma_k);
  }
  //foraging returns model  
  for ( i in 1:M ) {
         real m = log( alpha * phi[ID_i[i]] * psi[i]);
         returns[i] ~ lognormal( m , sigma );
  }
}
