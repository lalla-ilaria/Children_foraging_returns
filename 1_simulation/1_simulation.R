sim_data <- function (N , M , Q = 100,   #number of 
          zero = T, 
          b_a = 1,
          g_a = 1,
          g_k = 1,
          g_h = 1,
          g_g = 1,
          l = 1,
          alpha_success = 0.2,
          alpha_returns = 1,
          add_to_b = 0,
          max_age = age_plot){
  #parameters
  xi <-  rexp(1,l)
  sigma <- rexp(1, 4)
  
  #Simulate individual traits
  AGE <- runif(N,3,max_age)
  K <- abs( (AGE/mean(AGE)) + rnorm(N, 0, 0.3))#trait
  H <- add_to_b + abs( (AGE/mean(AGE)) + rnorm(N, 0, 0.3))#trait
  G <- abs( (AGE/mean(AGE)) + rnorm(N, 0, 0.3))#trait
  
  ########  
  ###Items
  ########
  # M items, each has unique difficulty (b) and discrimination (a)
  a <- rexp(Q, 1.5)
  b <- rnorm(Q, 0, 1)
  
  ##########
  ###Answers
  ##########
  #create a matrix to store answers
  Y <- matrix(NA,nrow=N,ncol=Q)
  #per each individual and question, find the probability p of a correct answer
  for ( i in 1:N ) for( j in 1:Q ) {
    p <- a[j]*( K[i] - b[j] )
    #assign correct or incorrect answer depending on prob of correct answer
    Y[i,j] <- rbern(1, inv_logit(p))
  }
  
  
  
  phi <- (1-exp( -b_a * AGE/mean(AGE) ) )^ g_a *
                  (K/mean(K)) ^ g_k *
                  (H/mean(H)) ^ g_h *
                  (G/mean(G)) ^ g_g

#simulate trip properties
  if (zero == T){
    #if zero are inflated, we're dealing with traps, hence time is in days
    L <- round ( runif(M, 1, 7) )#labor, duration of foraging trip (day)
  } else {
    #else, the trip is for shellfish collection, and the duration is in minutes
    L <- abs(rnorm(M, 150, 100)) #labor, duration of foraging trip (min)
  }
  #create matrices to save trip effects
  psi <- vector("numeric", length = M) 
  psi <-    L/mean(L) ^ xi
  ID_ind <- sample (1:N, size = M, replace = T)#assign trip to individual
  
  #calculate per datapoint
  lambda_pois <- vector("numeric", length = M) 
  S <- vector("numeric", length = M) 
  R <- vector("numeric", length = M) 
  #generates bernoulli distributed data - as if outcome data for traps were observations of traps
  # for(i in 1:M){
  #       if( zero == F) S <- rep(1, M) else {
  #           p[i] <- 1 - exp ( - abs(alpha_success) * phi[ID_ind[i]] * psi[i] )
  #           S[i] <- rbern(1, p[i])
  #           }     
  #       m <- log(alpha_returns * phi[ID_ind[i]] * psi[i])
  #       R[i] <- S[i] * rlnorm (1, m, 0.2)
  # }
  #generates poisson distributed data - as if outcome data were n of captures per trap (which is actual outcome)
  for(i in 1:M){
    if( zero == F) S <- rep(1, M) else {
      lambda_pois[i] <-  abs(alpha_success) * phi[ID_ind[i]] * psi[i] 
      S[i] <- rpois(1, lambda_pois[i])
    }     
    m <- log(alpha_returns * phi[ID_ind[i]] * psi[i])
    R[i] <- S[i] * rlnorm (1, m, 0.2)
  }
  has_knowledge <- rbern(N, 0.8)
  has_height <- rbern(N, 0.8)
  has_grip <- rbern(N, 0.8)

  return( list (N = N, #number of kids
                M = M, #number of trips
                Q = Q, #number of questions
                age = AGE, #ages
                sex = 1 + rbern(N, 0.5),
                knowledge = K, #knowledge
                height = H, #body
                grip = G, #body
                has_knowledge = has_knowledge,
                has_height = has_height,
                has_grip = has_grip,
                answers = Y, #answers to questionnaire
                duration = L, #duration  of trips
                ID_ind = ID_ind, #child per trip
                p = p, #prob success, for testing
                phi = phi, #for testing
                success = S,
                returns = R  #returns amount
                ) )
}





