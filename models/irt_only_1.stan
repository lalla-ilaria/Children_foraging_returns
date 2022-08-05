data {
//   int<lower=0> N;
//   int age_int[N];       //individual age
// 	int sex[N];       //individual sex
// 	int ID_k[N];      //index of forager/return for knowledge data
// 	int knowledge_impute[N]; //whether knowledge data are available or not
	//knowledge data
	int W;            //number of individuals in knowledge data
	int Q;            //number of questions
	int O;            //n ages prior drichlet
	int age_irt[W];   //age of individuals
	int sex_irt[W];   //sex of indiviudals
  int answers[W,Q]; //answers to freelist
	vector[O-1] prior_dirichlet; //prior drichlet
}

parameters {
	//knowledge parameters
	real<upper=0> omega;
	vector[W] iota_irt;       // individual intercepts on knowledge
  real<lower=0> sigma_irt;  //sd for iota_irt
  vector<lower=0>[2] ro_age; // coefficient relating age to knowledge, one per sex
  vector<lower=0>[Q] a;     //discrimination of questions
  vector[Q] b;     //difficulty of questions
  simplex[O-1] delta ;       //age specific effects
}
transformed parameters{
  //additional parameters
  vector[W] knowledge;
  //vector[N] knowledge_merged;
  vector[O] delta_j;
  delta_j  = append_row(1, delta);

  //knowledge estimation
  for ( i in 1:W )   knowledge[i] = omega + iota_irt[i] * sigma_irt +                    //individual interecepts -absorbs residual variation
                                    ro_age[sex_irt[i]] * sum (delta_j[ 1 : age_irt[i]] );//effect of age - sex specific
  //impute knowledge from irt if missing data
  // for(i in 1:N) {
  //   if (knowledge_impute[i] == 0) knowledge_merged[i] = knowledge[ID_k[i]]; 
  //   if (knowledge_impute[i] == 1) knowledge_merged[i] = 1;//omega + ro_age[sex[i]] * sum (delta_j[ 1 : age_int[i]] ); 
  // }//i
}
model {
	omega ~ normal( -5, 3)T[,0];
  for(i in 1:Q) a[i] ~ normal(0, 1) T[0,]; //value constrained above zero
	b ~ normal(0,2);
  //knowledge ~ normal(0,2);
	iota_irt ~ normal(0,1);
  for (i in 1:2) ro_age[i] ~ normal( 3 , 2 ) T[0,];
  delta ~ dirichlet( prior_dirichlet );
  sigma_irt ~ exponential(1);
  for(i in 1:W){
		vector[Q] p_irt = a .* (knowledge[i]-b);
		target += bernoulli_logit_lpmf(answers[i,]|p_irt);
	}
}

