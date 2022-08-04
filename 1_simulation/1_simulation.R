sim_data <- function (N , M , Q = 0,   #number of 
          zero = T, 
          b_a = 1,
          g_a = 1,
          g_k = 1,
          g_b = 1,
          l = 1,
          alpha = 0.3,
          add_to_b = 0){
  #parameters
  lambda <-  rexp(1,l)
  sigma <- rexp(1, 3)
  
  #Simulate individual traits
  AGE <- runif(N,3,20)
  K <- abs( (AGE/mean(AGE)) + rnorm(N, 0, 0.3))#trait
  B <- add_to_b + abs( (AGE/mean(AGE)) + rnorm(N, 0, 0.3))#trait
  
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
                  (B/mean(B)) ^ g_b

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
  psi <-    L/mean(L) ^ lambda
  ID_ind <- sample (1:N, size = M, replace = T)#assign trip to individual
  
  #calculate per datapoint
  p <- vector("numeric", length = M) 
  S <- vector("numeric", length = M) 
  R <- vector("numeric", length = M) 
  for(i in 1:M){
        if( zero == F) S <- rep(1, M) else {
            p[i] <- 1 - exp ( - abs(alpha) * phi[ID_ind[i]] * psi[i] )
            S[i] <- rbern(1, p[i])
            }     
        m <- log(alpha * phi[ID_ind[i]] * psi[i])
        R[i] <- S[i] * rlnorm (1, m, sigma)
  }

  return( list (N = N, #number of kids
                M = M, #number of trips
                Q = Q, #number of questions
                A = AGE, #ages
                K = K, #knowledge
                B = B, #body
                Y = Y, #answers to questionnaire
                L = L, #duration  of trips
                ID_ind = ID_ind, #child per trip
                p = p, #prob success, for testing
                phi = phi, #for testing
                S = S,
                R = R  #returns amount
                ) )
}





