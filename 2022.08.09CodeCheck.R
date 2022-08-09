library(rethinking)
library(dplyr)
library(rlist)
#load data
real_data <- list.load("2_data_preparation/processed_data.RData")
#define colors & stuff for plotting
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"
seq_trait <- seq(0,3,0.001)
age_plot <- 40

#############################
#PREPARE DATA
#############################
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells
d_shell_k <- real_data$shell_k

#remove people for whom we have only anthropometric data
d_shellppl <- d_shellppl %>% filter(!data == "anthropometrics")
d_shellppl <- d_shellppl %>% filter(!data == "knowledge")
d_shell_k <- d_shell_k[1:nrow(d_shellppl),]

#add index variables
#index and sort all individuals so we can loop across them
d_shellppl$index_id <- as.integer(as.factor(d_shellppl$anonymeID))
d_shellppl <- d_shellppl[order(d_shellppl$index_id),]
#add index for individuals in the foraging data
for ( i in 1:nrow(d_shells)){
  d_shells$index_id[i] <- d_shellppl$index_id[which ( d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
#sort knowledge data
d_shell_k <- d_shell_k[ order(row.names(d_shell_k)), ]

##############################
#FIT MODEL
##############################
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  has_foraging = ifelse(d_shellppl$data == "shells", 1, 0),
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #knowledge data
  Q = ncol(d_shell_k),                        #n items in freelist
  answers = d_shell_k                         #all answers from freelist
)
dat_shells[["answers"]][is.na(dat_shells[["answers"]])] <- -999

m_draft <- cstan( file= "models/3_draft.stan" , data=dat_shells , 
                      chains=3, cores = 3, iter = 500 )

#########################
#check outcomes
########################
#model fit
#(and does it badly)
precis(m_draft)
post <- extract.samples(m_draft)

#make plots
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)

#check standardized  knowledge
plot(NULL, xlim = c(1,dat_shells$N), ylim = c(0, 4), 
     xlab = "index", ylab = "standardized knowledge")
for (i in 1:dat_shells$N) {
  points ( rep(i, 750 ),
    post$knowledge_st[,i] ,
    col = rep(ifelse(dat_shells$has_knowledge[i] == 1, col.alpha("deepskyblue", 0.3), col.alpha("darkgoldenrod1", 0.3))))
}
points(1:dat_shells$N, apply(post$knowledge_st, 2, mean), 
     ylim = c(0, 4), pch = 19, col = presence_col)
for (i in 1:dat_shells$N) {
  lines ( rep(i, 2 ),
    PI(post$knowledge_st[,i]) ,
    col = rep(ifelse(dat_shells$has_knowledge[i] == 1, "deepskyblue4", "darkorange3")))
}

#check knowledge by age
plot(d_shellppl$knowledge, apply(post$knowledge_st, 2, mean))
year_eff <- apply(post$delta_j[,], 1, cumsum)
year_seq <- seq(0,4,0.2)
plot(x = dat_shells$age, 
     y = apply(post$knowledge, 2, mean), 
     xlim = c(0,4),
     xlab = "Age" , 
     ylab = "",
     cex.lab=1.8 , 
     cex.axis=1.8 ,
     pch = ifelse(dat_shells$has_knowledge == 1, 19, 1) , 
     cex = 1.5, 
     col =  alpha( sex_col , 0.6 )  )
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,1] * ( 1 - exp(-post$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:250) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,2] * ( 1 - exp(-post$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}

#plot of returns
post_s <- post
  phi <-  exp(mean(post_s$iota) ) * (
      (1-exp(- mean(post_s$beta) * seq_trait  )) ^ mean(post_s$gamma))
  psi <-   (mean(dat_shells$duration)) ^ mean(post_s$xi) * 
      exp(mean(dat_shells$tide) * mean(post_s$tau))
  R <- exp (  log(mean(post_s$alpha) * phi * psi) + 
                (mean(post_s$sigma)^2 /2))
  samp_data <- rlnorm(length(seq_trait),  
                      log(mean(post_s$alpha) * phi * psi), 
                      mean(post_s$sigma))
  
  plot(jitter(seq_trait) * mean(d_shellppl$age), samp_data, 
       xlim = c(0,age_plot), ylim = c(0, max(dat_shells$returns)+1), 
       xlab = "Age", ylab = "kg shellfish",
       pch = 16, col = col.alpha("orange", 0.2))
  for(i in 1:150){
    phi <-  exp(apply(post_s$iota,1,mean )[i] ) * (
      (1-exp(- post_s$beta[i] * seq_trait  )) ^ post_s$gamma[i]
    )
    psi <-   (mean(dat_shells$duration)) ^ post_s$xi[i] * 
      exp(mean(dat_shells$tide) * post_s$tau[i])
    R <- exp (  log(post_s$alpha[i] * phi * psi) + 
                  (post_s$sigma[i]^2 /2))
    lines( seq_trait * mean(d_shellppl$age),  R, 
           col = col.alpha(shellcol, 0.2), lwd = 1)
  }
  points(jitter(d_shellppl$age[d_shells$index_id], amount = 0.4), 
         as.numeric(d_shells$returns)/1000,
         pch = 16, cex = 0.8, col = col.alpha(othercol, 0.4))

