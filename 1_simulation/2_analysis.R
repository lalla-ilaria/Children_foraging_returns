
##################################################
#TEST FIT MODELS AGE ONLY
##################################################

#check that the model can fit to data
#note that the simulation generates return data from a lognormal distribution, so that the data can be generated in a variety of orders of magnitude and have to be re-scaled before being passed to the model
#SHELLS
d <- sim_data(50, 100, zero = F)#all parameters are 0.3
plot(1:d$M, d$returns)
order_of_magnitude <- ifelse(max(d$returns) <= 40, 0, floor(log(max(d$returns), 10)))
d$returns <- d$returns/10
#create data frame
dat_shells <- list(
  N = d$N,
  M = d$M,
  age = d$age/mean(d$age),
  returns = d$returns/ 10 ^ order_of_magnitude,
  duration = d$duration/mean(d$duration),
  tide = rep(0,d$M),
  ID_i= d$ID_ind
)
m_shell_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells , chains=3, cores = 3, iter = 1000 )

d <- sim_data(50, 100, zero = T)#all parameters are 0.3
dat_traps <- list(
  N = d$N,                       #n individuals in total sample
  M = d$M,                         #n trip/person
  ID_i= d$ID_ind,                    #index of person of trip 
  success = d$success,                 #whether trap captured something
  age = (d$age / mean(d$age)),
  duration = d$duration/mean(d$duration)
)
m_trap_age <- cstan( file= "models/1_trap_age_poisson.stan" , data=dat_traps , chains=3, cores = 3, iter = 1000 )

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

png("plots/validate_model/fit_to_data_validation.png", height = 2, width = 4, units = "in", res = 500)
par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
#shells
phi <-  mean(post_s$iota) +
  median(post_s$gamma) * log(1-exp(- mean(post_s$beta) * seq_trait  )) 
psi <-  median(post_s$xi) * (mean(log(dat_shells$duration))) + 
  median(post_s$tau)* mean(dat_shells$tide) 
samp_data <- rlnorm(length(seq_trait),  
                    mean(log(post_s$alpha)) + phi + psi, 
                    mean(post_s$sigma))


plot(jitter(seq_trait)* mean(d$age) , samp_data, 
     xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
     xlab = "Age", ylab = "kg shellfish",
     pch = 16, col = col.alpha("orange", 0.2))
for(i in 1:150){
  phi <-  apply(post_s$iota,1,mean )[i] +
    post_s$gamma[i] * log(1-exp(- post_s$beta[i] * seq_trait  )) 
  psi <-  post_s$xi[i] * (mean(log(dat_shells$duration))) +  
    post_s$tau[i]* mean(dat_shells$tide)
  R <- exp (  log(post_s$alpha)[i] + phi + psi + 
                (post_s$sigma[i]^2 /2))
  lines( seq_trait* mean(d$age) ,  R, 
         col = col.alpha(shellcol, 0.2), lwd = 1)
}

points(dat_shells$age[dat_shells$ID_i]* mean(d$age), dat_shells$returns,
       pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))
text(1, max(dat_shells$returns)+0.5, "A")

#traps 
phi <-  mean(post_t$iota) +
  mean(post_t$gamma) * log(1-exp(- mean(post_t$beta) * seq_trait  )) 
psi <-  mean(post_t$xi) * (mean(log (dat_traps$duration))) 
lambda <- mean(post_t$alpha) * exp(phi) * exp(psi)
samp_data <- rpois(length(seq_trait),  lambda)

#NB making plot with 13212 only (one of best hunters) to make plot with optimal situation to raise curve off zero
plot(jitter(seq_trait) * mean(d$age), samp_data, 
     xlab = "Age", ylab = "p capture",
     xlim = c(0,age_plot), ylim = c(0, 3.1), 
     pch = 16, col = col.alpha("lawngreen", 0.3))
#with average actor and average time 
for(i in 1:150){
  phi <-  apply(post_t$iota,1,mean )[i] +
    post_t$gamma[i] * log(1-exp(- post_t$beta[i] * seq_trait))
  psi <-  post_t$xi[i] * mean(log(dat_traps$duration))
  lambda <- post_t$alpha[i] * exp(phi) * exp(psi)
  lines( seq_trait * mean(d$age),  lambda,
         col = col.alpha(trapcol, 0.2), lwd = 1)
}
points(jitter(dat_traps$age[dat_traps$ID_i]* mean(d$age), amount = 0.5), 
       jitter(dat_traps$success, amount = 0.02), 
       pch = 16,  
       col = col.alpha(othercol, 0.6))
text(1, 0.95, "B")
dev.off()



##################################################
#FIT MODELS all data - check imputing missing data 
##################################################


d <- sim_data(50, 100, zero = F, g_g = 2, g_h = 2, g_k = 2)#all parameters are 0.3
#plot(1:d$M, d$returns)
order_of_magnitude <- ifelse(max(d$returns) <= 40, 0, floor(log(max(d$returns))))
#plot(d$age[d$ID_ind], d$returns/ 10 ^ order_of_magnitude)
dat_shells <- list(
  #foraging data
  N = d$N,                       #n individuals in total sample
  M = d$M,                         #n trip/person
  ID_i= d$ID_ind,                    #index of person of trip 
  returns = d$returns/ 10 ^ order_of_magnitude,#amount of shells in kg
  age = d$age / mean(d$age),
  sex = d$sex , #make vector of sexes 1 = male 2 = female
  duration = 1 + d$duration/mean(d$duration),
  tide = rep(0,d$M),
  #height data
  has_height = d$has_height,# #vector of 0/1 for whether height has to be imputed
  height = d$height/mean(d$height), 
  min_height = 0,#average height of newborn as intercept in height model
  #grip data
  has_grip = d$has_grip,# #vector of 0/1 for whether grip has to be imputed
  grip =  d$grip/mean(d$grip),
  #knowledge data
  has_knowledge = d$has_knowledge,# #vector of 0/1 for whether knowledge has to be imputed
  knowledge = d$knowledge/mean(d$knowledge),
  Q = d$Q,                        #n items in freelist
  answers = d$answers                         #all answers from freelist
)

m_shells_all <- cstan( file= "models/2_shells_all.stan" , data=dat_shells , 
                       chains=3, cores = 3, iter = 400 )


d <- sim_data(30, 500, zero = T, 
              g_g = 2, g_h = 2, g_k = 2,
              alpha_success = 0.02)#all parameters are 0.3

#TRAPS
dat_traps <- list(
  #foraging data
  N = d$N,                       #n individuals in total sample
  M = d$M,                         #n trip/person
  ID_i= d$ID_ind,                    #index of person of trip 
  success = d$success,                 #whether trap captured something
  age = d$age / mean(d$age),
  sex = d$sex , #make vector of sexes 1 = male 2 = female
  duration = d$duration/mean(d$duration),
  #height data
  has_height = d$has_height,# #vector of 0/1 for whether height has to be imputed
  height = d$height/mean(d$height),
  min_height = 0,#average height of newborn as intercept in height model
  #grip data
  has_grip = d$has_grip,# #vector of 0/1 for whether grip has to be imputed
  grip = d$grip/mean(d$grip),
  #knowledge data
  has_knowledge = d$has_knowledge,# #vector of 0/1 for whether knowledge has to be imputed
  knowledge = d$knowledge/mean(d$knowledge),
  Q = d$Q,                        #n items in freelist
  answers = d$answers                         #all answers from freelist
)



m_traps_all <- cstan( file= "models/2_traps_all_poisson.stan" , data=dat_traps , 
                      chains=3, cores = 3, iter = 400 )
post_s <- extract.samples(m_shells_all)
post_t <- extract.samples(m_traps_all)

###################
#missing data
###################
#make plots

#check standardized  knowledge
png("plots/validate_model/knowledge_estimates_validation.png", height = 3, width = 4, units = "in", res = 500)
par(mfrow = c(1,2),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
#shells
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", shellcol, "deepskyblue4")
plot(NULL, xlim = c(0,2.5), ylim = c(0, 3.5), 
     xlab = "simulated knowledge", ylab = "irt knowledge")
for (i in 1:dat_shells$N) {
  points ( rep(dat_shells$knowledge[i], nrow(post_s$knowledge) ),
           post_s$knowledge[,i] , pch = 19, 
           cex = 0.3,
           col = rep(ifelse(dat_shells$has_knowledge[i] == 1, col.alpha("orange", 0.1), col.alpha("deepskyblue", 0.1))))
}
points(dat_shells$knowledge, apply(post_s$knowledge, 2, mean), 
       ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_shells$N) {
  lines ( rep(dat_shells$knowledge[i], 2 ),
          PI(post_s$knowledge[,i]) ,
          col = rep(ifelse(dat_shells$has_knowledge[i] == 1, shellcol, "deepskyblue4")))
}
#traps
sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_traps$has_knowledge == "1", trapcol, "deepskyblue4")
plot(NULL, xlim = c(0,2.5), ylim = c(0, 3.5), 
     xlab = "simulated knowledge", ylab = "irt knowledge")
for (i in 1:dat_traps$N) {
  points ( rep(dat_traps$knowledge[i], nrow(post_t$knowledge) ),
           post_t$knowledge[,i] , pch = 19, 
           cex = 0.3,
           col = rep(ifelse(dat_traps$has_knowledge[i] == 1, col.alpha("lawngreen", 0.1), col.alpha("deepskyblue", 0.1))))
}
points(dat_traps$knowledge, apply(post_t$knowledge, 2, mean), 
       ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_traps$N) {
  lines ( rep(dat_traps$knowledge[i], 2 ),
          PI(post_t$knowledge[,i]) ,
          col = rep(ifelse(dat_traps$has_knowledge[i] == 1, trapcol, "deepskyblue4")))
}

dev.off()

#check traits by age shells
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", shellcol, "deepskyblue4")
png("plots/validate_model/missing_data_validation_shells.png", height = 5, width = 8, units = "in", res = 500)
par(mfrow = c(1,3),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
year_seq <- seq(0,4,0.2)
plot(x = dat_shells$age, 
     y = apply(post_s$height_merged, 2, mean), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "Height",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_height == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq ,  
        y = dat_shells$min_height + post_s$kappa[i,1] * ( 1 - exp(-post_s$chi[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq ,  
        y = dat_shells$min_height + post_s$kappa[i,2] * ( 1 - exp(-post_s$chi[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_shells$age, 
     y = apply(post_s$grip_merged, 2, mean), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "grip",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_grip == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_s$epsilon[i,1] * ( 1 - exp(-post_s$upsilon[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_s$epsilon[i,2] * ( 1 - exp(-post_s$upsilon[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_shells$age, 
     y = apply(post_s$knowledge, 2, median), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_s$omega[i] + post_s$omicron[i,1] * ( 1 - exp(-post_s$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_s$omega[i] + post_s$omicron[i,2] * ( 1 - exp(-post_s$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
dev.off()



#make plots
sex_col <- ifelse(dat_traps$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_traps$has_knowledge == "1", "deepskyblue4", "darkorange3")

png("plots/validate_model/missing_data_validation_traps.png", height = 5, width = 8, units = "in", res = 500)
par(mfrow = c(1,3),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
year_seq <- seq(0,4,0.2)
plot(x = dat_traps$age, 
     y = apply(post_t$height_merged, 2, mean), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "Height",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_height == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq ,  
        y = dat_traps$min_height + post_t$kappa[i,1] * ( 1 - exp(-post_t$chi[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq ,  
        y = dat_traps$min_height + post_t$kappa[i,2] * ( 1 - exp(-post_t$chi[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_traps$age, 
     y = apply(post_t$grip_merged, 2, mean), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "grip",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_grip == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_t$epsilon[i,1] * ( 1 - exp(-post_t$upsilon[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq ,  
        y = post_t$epsilon[i,2] * ( 1 - exp(-post_t$upsilon[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

plot(x = dat_traps$age,
     y = apply(post_t$knowledge, 2, median), 
     xlim = c(0,2),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_traps$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:100) {
  lines(x = year_seq,
        y = post_t$omega[i] + post_t$omicron[i,1] * ( 1 - exp(-post_t$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:100) {
  lines(x = year_seq, 
        y = post_t$omega[i] + post_t$omicron[i,2] * ( 1 - exp(-post_t$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
dev.off()


##############################################################
#FIT MODELS all data - model validation and parameter recovery 
##############################################################

#define parameter sweep
pars <- data.frame( g_g = c(2, 1, 3,-1, 4,-3, -2, 0),
                    g_h = c(2, 4,-4,-3, 0, 1, 3, -1),
                    g_k = c(2,-2, 1, 0,-3, 4, -1, 3))

out_shells <- list()
out_traps <- list()
phi_trait <- list()
shell_posteriors <- list()
traps_posteriors <- list()
rates_traps <- c(0.01, 0.1, 0.5, 0.8) #vector(mode = "numeric",length = nrow(pars))
rates_shells <- c(0.1, 0.5, 1, 2)
alphas <-vector(mode = "numeric",length = nrow(pars))

for(a in 1:4){
  png(paste("plots/validate_model/validation_parameters_recovery", rates_traps[a], ".png", sep = ""), height = 5, width = 8, units = "in", res = 500)
  par(mfrow = c(1,1),mgp = c(1.5, 0.5, 0), mar = c(2.5, 2.5, 2, 1) + 0.1)
  plot(NULL, 
       xlim = c( min(pars) - 1, max(pars) +1 ),
       ylim = c( min(pars) - 2, max(pars) +2 ), 
       xlab = "'real' parameter", ylab = "recovered parameter")
  legend(min(pars)-0.8, max(pars)+1.5, legend = c(expression(~theta), expression(~eta),expression(~zeta)), col = othercol, pch=15:17)
  
  for ( j in 1:nrow(pars)){
    #SHELLS
    d <- sim_data(50, 100, zero = F, g_g = pars[j,"g_g"], g_h = pars[j,"g_h"], g_k = pars[j,"g_k"], alpha_returns = rates_shells[a])#all parameters are 0.3
    order_of_magnitude <- ifelse(max(d$returns) <= 40, 0, floor(log(max(d$returns))))
    dat_shells <- list(
      #foraging data
      N = d$N,                       #n individuals in total sample
      M = d$M,                         #n trip/person
      ID_i= d$ID_ind,                    #index of person of trip 
      returns = d$returns/ 10 ^ order_of_magnitude,#amount of shells in kg
      age = d$age / mean(d$age),
      sex = d$sex , #make vector of sexes 1 = male 2 = female
      duration = 1 + d$duration/mean(d$duration),
      tide = rep(0,d$M),
      #height data
      has_height = d$has_height,# #vector of 0/1 for whether height has to be imputed
      height = d$height/mean(d$height),
      min_height = 0,#average height of newborn as intercept in height model
      #grip data
      has_grip = d$has_grip,# #vector of 0/1 for whether grip has to be imputed
      grip = d$grip/mean(d$grip),
      #knowledge data
      has_knowledge = d$has_knowledge,# #vector of 0/1 for whether knowledge has to be imputed
      #knowledge = d$knowledge/mean(d$knowledge)
      Q = d$Q,                        #n items in freelist
      answers = d$answers                         #all answers from freelist
    )
    
    m_shells_all <- cstan( file= "models/2_shells_all.stan" , data=dat_shells , 
                           chains=3, cores = 3, iter = 400 )
    
    
    d <- sim_data(30, 1000, zero = T, 
                  g_g = pars[j,"g_g"], g_h = pars[j,"g_h"], g_k = pars[j,"g_k"],
                  alpha_success = rates_traps[a])
    
    #TRAPS
    dat_traps <- list(
      #foraging data
      N = d$N,                       #n individuals in total sample
      M = d$M,                         #n trip/person
      ID_i= d$ID_ind,                    #index of person of trip 
      success = d$success,                 #whether trap captured something
      age = d$age / mean(d$age),
      sex = d$sex , #make vector of sexes 1 = male 2 = female
      duration = d$duration/mean(d$duration),
      #height data
      has_height = d$has_height,# #vector of 0/1 for whether height has to be imputed
      height = d$height/mean(d$height),
      min_height = 0,#average height of newborn as intercept in height model
      #grip data
      has_grip = d$has_grip,# #vector of 0/1 for whether grip has to be imputed
      grip = d$grip/mean(d$grip),
      #knowledge data
      has_knowledge = d$has_knowledge,# #vector of 0/1 for whether knowledge has to be imputed
      #knowledge = d$knowledge/mean(d$knowledge)
      Q = d$Q,                        #n items in freelist
      answers = d$answers                         #all answers from freelist
    )
    
    
    
    m_traps_all <- cstan( file= "models/2_traps_all_poisson.stan" , data=dat_traps , 
                          chains=3, cores = 3, iter = 400 )
    
    post_s <- extract.samples(m_shells_all)
    post_t <- extract.samples(m_traps_all)
    shell_posteriors[[j]] <- post_s
    traps_posteriors[[j]] <- post_t
    points(pars[j,1], mean(post_s$theta_g), col = col.alpha(shellcol, 0.6), pch = 15)
    points(pars[j,2], mean(post_s$eta_h), col = col.alpha(shellcol, 0.6), pch = 16)
    points(pars[j,3], mean(post_s$zeta_k), col = col.alpha(shellcol, 0.6), pch = 17)
    points(pars[j,1], mean(post_t$theta_g), col = col.alpha(trapcol, 0.6), pch = 15)
    points(pars[j,2], mean(post_t$eta_h), col = col.alpha(trapcol, 0.6), pch = 16)
    points(pars[j,3], mean(post_t$zeta_k), col = col.alpha(trapcol, 0.6), pch = 17)
    lines(rep(pars[j,1], 2), PI(post_s$theta_g), col = col.alpha(shellcol, 0.8))
    lines(rep(pars[j,2], 2), PI(post_s$eta_h), col = col.alpha(shellcol, 0.8))
    lines(rep(pars[j,3], 2), PI(post_s$zeta_k), col = col.alpha(shellcol, 0.8))
    lines(rep(pars[j,1], 2), PI(post_t$theta_g), col = col.alpha(trapcol, 0.8))
    lines(rep(pars[j,2], 2), PI(post_t$eta_h), col = col.alpha(trapcol, 0.8))
    lines(rep(pars[j,3], 2), PI(post_t$zeta_k), col = col.alpha(trapcol, 0.8))
    
    #CALCULATE DIFFERENCE BETWEEN ONE STANDARD DEVIATION OF SIMULATED DATA
    diffs_out <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 9))
    colnames(diffs_out) <- c("sg","sh","sk","sd","st","tg","th","tk","td")
    diffs_phi <- data.frame(matrix(nrow=length(post_s$alpha), ncol = 6))
    colnames(diffs_phi) <- c("sg","sh","sk","tg","th","tk")
    
    sds <- c( sd(log(dat_shells$grip)),
              sd(log(dat_shells$height)), 
              sd(apply(post_s$knowledge, 2, mean)), #sd(log(dat_shells$knowledge))
              sd(log(dat_shells$duration)),
              sd(dat_shells$tide))
    
    for (i in 1:5) {
      phi_min <-  apply(post_s$iota,1,mean )  +
        post_s$gamma * log(1-exp(- post_s$beta * 1.6  )) +
        post_s$theta_g*(mean(log(dat_shells$grip), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
        post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
        #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) - ifelse(i == 1, sds[i], 0))
        post_s$zeta_k*(mean(apply(post_s$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
      psi_min <-  post_s$xi * (mean(log(dat_shells$duration)) - ifelse(i == 4, sds[i], 0))+
        post_s$tau* (mean(dat_shells$tide) - ifelse(i == 5, sds[i], 0))
      min_k <- exp (log(post_s$alpha) + phi_min + psi_min +
                      (post_s$sigma^2 /2))
      phi_max <-  apply(post_s$iota,1,mean )  +
        post_s$gamma * log(1-exp(- post_s$beta * 1.6  )) +
        post_s$theta_g* (mean(log(dat_shells$grip), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
        post_s$eta_h* (mean(log(dat_shells$height), na.rm = TRUE) + ifelse(i == 2, sds[i], 0)) +
        #post_s$zeta_k* (mean(log(dat_shells$knowledge), na.rm = TRUE) + ifelse(i == 1, sds[i], 0))
        post_s$zeta_k* (mean(apply(post_s$knowledge, 2, mean)) + ifelse(i == 3, sds[i], 0))
      psi_max <-  post_s$xi * (mean(log(dat_shells$duration)) + ifelse(i == 4, sds[i], 0))+
        post_s$tau* (mean(dat_shells$tide) + ifelse(i == 5, sds[i], 0))
      max_k <- exp (log(post_s$alpha) + phi_max + psi_max +
                      (post_s$sigma^2 /2))
      diffs_out[i] <- max_k - min_k
      if( i %in% 1:3) diffs_phi[i] <- phi_max - phi_min
    }
    
    diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                      rep("grip", nrow(diffs_out)),
                                      rep("height", nrow(diffs_out)),
                                      rep("duration", nrow(diffs_out)),
                                      rep("tide", nrow(diffs_out))),
                        kg_shellfish =  c(diffs_out$sk, diffs_out$sg, diffs_out$sh, diffs_out$sd, diffs_out$st)
    )
    
    diffs_shells <- diffs %>% 
      mutate( variable = fct_reorder(.f = variable, .x = kg_shellfish, .fun = mean))
    
    
    sds <- c( sd(log(dat_traps$grip)),
              sd(log(dat_traps$height)), 
              sd(apply(post_t$knowledge, 2, mean)), #sd(log(dat_traps$knowledge))
              sd(log(dat_traps$duration)))
    
    for (i in 1:4) {
      phi_min <-  apply(post_t$iota,1,mean )  +
        post_t$gamma * log(1-exp(- post_t$beta * 1.6  )) +
        post_t$theta_g*(mean(log(dat_traps$grip), na.rm = TRUE) - ifelse(i == 1, sds[i], 0)) +
        post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) - ifelse(i == 2, sds[i], 0)) +
        #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) - ifelse(i == 3, sds[i], 0))
        post_t$zeta_k*(mean(apply(post_t$knowledge, 2, mean)) - ifelse(i == 3, sds[i], 0))
      psi_min <-  post_t$xi * (mean(log(dat_traps$duration)) - ifelse(i == 4, sds[i], 0))
      min_k <- 1 - exp ( - post_t$alpha * exp(phi_min) * exp(psi_min))
      phi_max <-  apply(post_t$iota,1,mean )  +
        post_t$gamma * log(1-exp(- post_t$beta * 1.6  )) +
        post_t$theta_g* (mean(log(dat_traps$grip), na.rm = TRUE) + ifelse(i == 1, sds[i], 0)) +
        post_t$eta_h* (mean(log(dat_traps$height), na.rm = TRUE) + ifelse(i == 2, sds[i], 0)) +
        #post_s$zeta_k* (mean(log(dat_traps$knowledge), na.rm = TRUE) ifelse(i == 1, sds[i], 0)) +
        post_t$zeta_k* (mean(apply(post_t$knowledge, 2, mean)) + ifelse(i == 3, sds[i], 0)) 
      psi_max <-  post_t$xi * (mean(log(dat_traps$duration)) + ifelse(i == 4, sds[i], 0))
      max_k <- 1 - exp ( - post_t$alpha * exp(phi_max) * exp(psi_max))
      diffs_out[i + 5] <- max_k - min_k
      if( i %in% 1:3) diffs_phi[i + 3] <- phi_max - phi_min
    }
    
    diffs <- data.frame(variable = c( rep("knowledge", nrow(diffs_out)),
                                      rep("grip", nrow(diffs_out)),
                                      rep("height", nrow(diffs_out)),
                                      rep("duration", nrow(diffs_out))),
                        p_success =  c(diffs_out$tk, diffs_out$tg, diffs_out$th, diffs_out$td)
    )
    
    diffs_traps <- diffs %>% 
      mutate( variable = fct_reorder(.f = variable, .x = p_success, .fun = mean))
    
    out_shells[[j]] <- ggplot(diffs_shells, aes(x = variable, y = kg_shellfish)) + 
      geom_hline(yintercept = 0, color = "grey90") +
      geom_violin(fill = col.alpha(shellcol, 0.6),
                  colour = shellcol) +
      labs(y = "kg shellfish difference", x = "")+
      ylim(-4, 6)+
      theme_classic()
    
    out_traps[[j]] <- ggplot(diffs_traps, aes(x = variable, y = p_success)) + 
      geom_hline(yintercept = 0, color = "grey90") +
      geom_violin(fill = col.alpha(trapcol, 0.6),
                  colour = trapcol) +
      labs(y = "p_success difference", x = "")+
      ylim(-1, 1)+
      theme_classic()
    
    v_diffs_phi <- data.frame(variable = c( rep("knowledge", nrow(diffs_phi)),
                                            rep("grip", nrow(diffs_phi)),
                                            rep("height", nrow(diffs_phi)),
                                            rep("knowledge", nrow(diffs_phi)),
                                            rep("grip", nrow(diffs_phi)),
                                            rep("height", nrow(diffs_phi))),
                              p_success =  c(diffs_phi$sk, diffs_phi$sg, diffs_phi$sh, diffs_phi$tk, diffs_phi$tg, diffs_phi$th),
                              foraging_type = c(rep("shell", 3*nrow(diffs_phi)), rep("trap", 3*nrow(diffs_phi))))
    
    
    phi_trait[[j]] <- ggplot(v_diffs_phi, aes(x = variable, y = p_success, fill = foraging_type , color = foraging_type) )+ 
      geom_hline(yintercept = 0, color = "grey90") +
      geom_violin(trim=FALSE)+
      scale_fill_manual("foraging_type",  limits= c("shell", "trap"), values = c( col.alpha(shellcol, 0.6), col.alpha(trapcol, 0.6 )),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
      scale_color_manual("foraging_type",  limits= c("shell", "trap"), values = c(shellcol, trapcol),guide = guide_legend())+ #gives colors as defined. Add "inland", and "#FFFFFF", for marine resorurces
      labs(x = "", y = "phi difference")+
      ylim(-10, 10)+
      theme_classic()+
      theme(legend.position="top")
    alphas[j] <- sum(d$success)
  }
  dev.off()
  
}

# for(j in 1:nrow(pars)){
#   plots <- align_plots(phi_trait[[j]], out_shells[[j]], align = 'v', axis = 'l')
#   # then build the bottom row
#   bottom_row <- plot_grid(plots[[2]], out_traps[[j]], labels = c('B', 'C'), label_size = 12)
#   
#   # then combine with the top row for final plot
#   plot_grid(plots[[1]], bottom_row, labels = c(paste(colnames(pars), pars[j,],collapse=","), ''), label_size = 12, ncol = 1)
#   ggsave(filename = paste("plots/validate_model/look_outcomes", j ,".png", sep = ""),
#          height = 14, width = 16, units = "cm", dpi = 500)
#   
# }
# 
