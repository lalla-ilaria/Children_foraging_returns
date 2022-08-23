data{
	//individual level data
	int N;                  //total number of individuals
	array[N] real age;      //individual age
	array[N] int sex;       //individual sex
	array[N] real height;   //individual height
	real min_height;        //height of newborn as intercept
	array[N] real grip;     //individual grip

	//data availability
	array[N] int has_knowledge;
	array[N] int has_height;
	array[N] int has_grip;

	// foraging trip data
	int M;                  //number of trip
	array[M] int ID_i;      //index of forager/return
	array[M] int success;  //success
	array[M] real duration; //length of trip

	//knowledge data
	int Q;                  //number of questions
	array[N,Q] int answers; //answers to freelist
}//data
	
parameters{
  //foraging parameters
  vector [N] iota;       //individual level random effect
  real<lower=0> sigma_i; //sd of iota
  real<lower=0> alpha;   //scaling intercept
  real<lower=0> beta;    //age effect
  real<lower=0> gamma;   //age elasticity
  real<lower=0> zeta_k;  //knowledge elasticity
  real<lower=0> eta_h;   //height elasticity
  real<lower=0> theta_g; //grip elasticity
  real xi;               //exponent for length trip

	//knowledge parameters
	real omega;
	vector[N] iota_knowledge;       //individual intercepts on knowledge
  real<lower=0> sigma_knowledge;  //sd for iota_knowledge
  vector<lower=0>[2] chi;   //max effect of age on knowledge, sex specific
  vector<lower=0>[2] ro_age;//sex specific age speed
  vector<lower=0>[Q] a;     //discrimination of questions
  vector<lower=0>[Q] b;     //difficulty of questions
  
  //height parameters
  real<lower=0> sigma_height; //sd for height estimation
  vector<lower=0>[2] kappa;   //max effect of age on height, sex specific
  vector<lower=0>[2] lambda;  //sex specific height speed
  vector<lower=0>[N] height_unobs;//vector of distributions for height - used when unobserved
  
  //grip parameters
  real<lower=0> sigma_grip;    //sd for height estimation
  vector<lower=0>[2] epsilon;  //max effect of age on grip, sex specific
  vector<lower=0>[2] upsilon;  //sex specific grip speed
  vector<lower=0>[N] grip_unobs;//vector of distributions for grip - used when unobserved
  real min_grip;
}//parameters

transformed parameters{
  //additional parameters
  vector[N] knowledge;
  vector[N] height_merged;
  vector[N] grip_merged;
  vector[N] phi;              //individual total effect
  vector[M] psi;              //trip effect
  
  //estimate knowledge for all individuals
  for ( i in 1:N ){
        knowledge[i] = omega +
                 iota_knowledge[i]*sigma_knowledge +
                 chi[sex[i]] * ( 1 - exp(-ro_age[sex[i]] * age[i]) );
  }

  //estimate height for all individuals
  for (i in 1:N){
    if(has_height[i] == 1) height_merged[i] = height[i];
    if(has_height[i] == 0) height_merged[i] = height_unobs[i];
  }
  
  //estimate grip for all individuals
  for (i in 1:N){
    if(has_grip[i] == 1) grip_merged[i] = grip[i];
    if(has_grip[i] == 0) grip_merged[i] = grip_unobs[i];
  }

  //estimate individual level parameters
  for(i in 1:N) {
    phi[i] = iota[i] * sigma_i + 
             gamma*log(1-exp(-beta * age[i])) + 
             eta_h*height_merged[i] + 
             theta_g*grip_merged[i] +
             zeta_k*knowledge[i];
    }
    
  //estimate trip level parameters  
  for(i in 1:M) {
    psi[i] =  xi*log(duration[i]); 
    }
}//transformed parameters

model{
  //priors for foraging parameters
  iota ~ normal(0,1);
  sigma_i ~ exponential(1);
  alpha ~ normal(0,1)T[0,];
  beta ~ exponential(1);
  gamma~ exponential(1);
  zeta_k~ exponential(1);
  eta_h~ exponential(1);
  xi ~ normal(0, 1);

  //priors for knowledge parameters
  omega ~ normal (-5,2);
  a ~ exponential(1);
	for(i in 1:Q) b[i] ~ normal(1,0.5);
	iota_knowledge ~ normal(0,1);
  chi ~ exponential(0.2);
  ro_age ~ exponential(1);
  sigma_knowledge ~ exponential(1);
  
  //priors for height parameters
  sigma_height ~ exponential(1);
  for (i in 1:N) height_unobs[i] ~ normal(1, 0.5)T[0,];
  kappa ~ exponential(1);
  lambda ~ exponential(1);
  
  //priors for grip parameters
  sigma_grip ~ exponential(1);
  for (i in 1:N) grip_unobs[i] ~ normal(1, 0.5)T[0,];
  epsilon ~ exponential(1);
  upsilon ~ exponential(1);
  min_grip ~ normal(0, 0.5);

	//irt across knowledge people to estimage age and sex specific knowledge
	for ( i in 1:N ) {
    if (has_knowledge[i] == 1){
		  for (j in 1:Q ) {
		    real p_irt = inv_logit(a[j]*(knowledge[i]-b[j]));
			  answers[i,j] ~ bernoulli( p_irt );
		  }
		}
	}
  
  //height model
  for (i in 1:N){
      real mu_height = min_height + kappa[sex[i]] * ( 1 - exp(-lambda[sex[i]] * age[i]) );
      height_merged[i] ~ normal(mu_height, sigma_height);

  }

  //grip model
  for (i in 1:N){
      real mu_grip = min_grip + epsilon[sex[i]] * ( 1 - exp(-upsilon[sex[i]] * age[i]) );
      grip_merged[i] ~ normal(mu_grip, sigma_grip);

  }
  
  //foraging returns model  
  for ( i in 1:M ) {
         real p = 1 - exp (- alpha * exp(phi[ID_i[i]]) * exp(psi[i]));
         success[i] ~ bernoulli(p); 
      }
}
