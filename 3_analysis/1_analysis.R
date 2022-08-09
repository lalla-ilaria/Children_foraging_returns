library(rethinking)
library(dplyr)
library(rlist)
real_data <- list.load("2_data_preparation/processed_data.RData")
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"

#d_knowledge <- list.load("../data/processed_data.RData")


##########################################################################
#PREPARE DATA -  AGE
##########################################################################
#donow4 scale variables so they are above 1 (one at a time, check what happens)
#SHELLS
#keep people for whom we have ages
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$age),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]
#add index variables
dc_shellppl$index_id <- as.integer(as.factor(dc_shellppl$anonymeID))
dc_shellppl <- dc_shellppl[order(dc_shellppl$index_id),]
dc_shells$index_id <- as.integer(as.factor(dc_shells$anonymeID))
dc_shell_k <- dc_shell_k[ order(as.factor(row.names(dc_shell_k))), ]
#create data frame
dat_shells <- list(
  N = nrow(dc_shellppl),
  M = nrow(dc_shells),
  Q = ncol(dc_shell_k),
  A = dc_shellppl$age / mean(dc_shellppl$age),
  R = as.numeric(dc_shells$returns)/1000,
  L = dc_shells$lenght_min/mean(dc_shells$lenght_min),
  tide = dc_shells$tide_avg_depth,
  ID_i= dc_shells$index_id
)

#TRAPS
#keep complete cases only
dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$age),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]
#add index variables
dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]
#create data frame
dat_traps <- list(
  N = nrow(dc_trapppl),
  M = nrow(dc_traps),
  Q = ncol(dc_trap_k),
  A = dc_trapppl$age / mean(dc_trapppl$age),
  S = as.numeric(dc_traps$success),
  L = dc_traps$lenght_hour/mean(dc_traps$lenght_hour),
  ID_i= as.integer(as.factor(as.character(dc_traps$anonymeID)))
)

#model with age only
m_shell_age <- cstan( file= "models/1_shell_age.stan" , data=dat_shells , chains=3, cores = 3 )
m_trap_age <- cstan( file= "models/1_trap_age.stan" , data=dat_traps , chains=3, cores = 3 )

##########################################################################
#PREPARE DATA - ALL
##########################################################################
#donow4 scale variables so they are above 1 (one at a time, check what happens)
#SHELLS
#keep complete cases only
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$knowledge),]
dc_shellppl <- dc_shellppl[complete.cases(dc_shellppl$height),]
dc_shellppl <- dc_shellppl[complete.cases(dc_shellppl$grip),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]
#add index variables
dc_shellppl$index_id <- as.integer(as.factor(dc_shellppl$anonymeID))
dc_shellppl <- dc_shellppl[order(dc_shellppl$index_id),]
dc_shells$index_id <- as.integer(as.factor(dc_shells$anonymeID))
dc_shell_k <- dc_shell_k[ order(as.factor(row.names(dc_shell_k))), ]
#create data frame
dat_shells <- list(
  N = nrow(dc_shellppl),
  M = nrow(dc_shells),
  Q = ncol(dc_shell_k),
  A = dc_shellppl$age / mean(dc_shellppl$age),
  Y = dc_shell_k,
  K = dc_shellppl$knowledge / mean(dc_shellppl$knowledge),
  H = dc_shellppl$height / mean(dc_shellppl$height),
  G = dc_shellppl$grip / mean(dc_shellppl$grip),
  R = as.numeric(dc_shells$returns)/1000,
  L = dc_shells$lenght_min/mean(dc_shells$lenght_min),
  tide = dc_shells$tide_avg_depth,
  ID_i= dc_shells$index_id
)

#TRAPS
#keep complete cases only
dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$knowledge),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]
#add index variables
dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]
#create data frame
dat_traps <- list(
  N = nrow(dc_trapppl),
  M = nrow(dc_traps),
  Q = ncol(dc_trap_k),
  A = dc_trapppl$age / mean(dc_trapppl$age),
  Y = dc_trap_k,
  K = dc_trapppl$knowledge / mean(dc_trapppl$knowledge),
  H = dc_trapppl$height / mean(dc_trapppl$height),
  G = dc_trapppl$grip / mean(dc_trapppl$grip),
  S = as.numeric(dc_traps$success),
  L = dc_traps$lenght_hour/mean(dc_traps$lenght_hour),
  ID_i= as.integer(as.factor(as.character(dc_traps$anonymeID)))
)

#model with all predictors
m_shell_all <- cstan( file= "models/2_shell_all.stan" , data=dat_shells , chains=3, cores = 3 )
m_trap_all <- cstan( file= "models/2_trap_all.stan" , data=dat_traps , chains=3, cores = 3 )


#CHECKS
#always check that knowledge is not hitting zero (maybe with a function?)

##########################################################################
#PREPARE DATA -  input knowledge - draft
##########################################################################
#SHELLS
#create data frame
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells
d_shell_k <- real_data$shell_k
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps
d_trap_k <- real_data$trap_k


#remove people for whom we have only anthropometric data
d_shellppl <- d_shellppl %>% filter(!data == "anthropometrics")
d_shellppl <- d_shellppl %>% filter(!data == "knowledge")
d_shell_k <- d_shell_k[1:nrow(d_shellppl),]
d_trapppl <- d_trapppl %>% filter(!data == "anthropometrics")
d_trapppl <- d_trapppl %>% filter(!data == "knowledge")
d_trap_k <- d_trap_k[1:nrow(d_trapppl),]



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

#same for traps
d_trapppl$index_id <- as.integer(as.factor(d_trapppl$anonymeID))
d_trapppl <- d_trapppl[order(d_trapppl$index_id),]
for ( i in 1:nrow(d_traps)){
  d_traps$index_id[i] <- d_trapppl$index_id[which ( d_trapppl$anonymeID == d_traps$anonymeID[i])]
}
d_trap_k <- d_trap_k[ order(row.names(d_trap_k)), ]

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  has_foraging = ifelse(d_shellppl$data == "shells", 1, 0),#[1:40],
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),#[1:40], #vector of 0/1 for whether knowledge has to be imputed
  returns = as.numeric(d_shells$returns)/1000,
  age = (d_shellppl$age / mean(d_shellppl$age)),#[1:40],
  age_int = d_shellppl$age,#[1:40],
  sex = ifelse(d_shellppl$sex == "m", 1, 2),#[1:40], #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #knowledge data
  Q = 100,#ncol(d_shell_k),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  answers = d_shell_k[, 1:100], #all answers from freelist
  prior_dirichlet = rep( 0.5, length (1:56 ) -1 ) #prior for Dirichlet
)
for (i in 1:length(dat_shells)) dat_shells[[i]][is.na(dat_shells[[i]])] <- -999

#try tightening priors
#making a positive version of knowledge
#connect to foraging stuff
m_shell_irt <- cstan( file= "models/2_shell_k3.stan" , data=dat_shells , 
                      chains=3, cores = 3, iter = 500 )

#m <- cmdstan_model("models/2_shell_k2.stan" )
#standardizing knowledge IN the models???
post <- extract.samples(m_shell_irt)
sex_col <- ifelse(dat_shells$sex == "1", boycol, girlcol)
presence_col <- ifelse(dat_shells$has_knowledge == "1", "deepskyblue4", "darkorange3")
plot(NULL, xlim = c(1,dat_shells$N), ylim = c(0, 4))
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
for (i in 1:150) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,1] * ( 1 - exp(-post$ro_age[i,1] * year_seq)) ,  
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:150) {
  lines(x = year_seq ,  
        y = post$omega[i] + post$chi[i,2] * ( 1 - exp(-post$ro_age[i,2] * year_seq)) ,  
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}


for (i in 1:150) {
  lines(x = 1:nrow(year_eff),  
        y = post$omega[i] + post$ro_age[i,1] * year_eff[,i], 
        type = "l", 
        col = col.alpha( boycol, alpha = 0.1))}
for (i in 1:150) {
  lines(x = 1:nrow(year_eff),  
        y = post$omega[i] + post$ro_age[i,2] * year_eff[,i], 
        type = "l", 
        col = col.alpha( girlcol, alpha = 0.1))}
lines( x = 1:nrow(year_eff),  
       y = mean(post$omega) + mean(post$ro_age[,1]) * apply(year_eff, 1, mean), 
       type = "l", 
       col = col.alpha( boycol, alpha = 1), lwd = 1.5 )
lines( x = 1:nrow(year_eff),  
       y = mean(post$omega) + mean(post$ro_age[,2]) * apply(year_eff, 1, mean), 
       type = "l", 
       col = col.alpha( girlcol, alpha = 1), lwd = 1.5 )

#m_shell_all <- cstan( file= "models/2_shell_knowledgeinput.stan" , data=dat_shells , chains=3, cores = 3 )

#TRAPS
#keep complete cases only
d_trapppl <- real_data$trap_ppl
d_traps <- real_data$traps
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]
#add index variables
dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]
#create data frame
dat_traps <- list(
  N = nrow(dc_trapppl),
  M = nrow(dc_traps),
  Q = ncol(dc_trap_k),
  A = dc_trapppl$age / mean(dc_trapppl$age),
  S = as.numeric(dc_traps$success),
  L = dc_traps$lenght_hour/mean(dc_traps$lenght_hour),
  ID_i= as.integer(as.factor(as.character(dc_traps$anonymeID)))
)

dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = d_shell_k_ppl$index_id, #make vector of length W where is reported index ID of individuals in foraging df
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = ncol(d_knowledge$Y_l),
  O = 57 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_knowledge$Y_l, #all answers from freelist
  prior_dirichlet = rep( 0.5, length (0:40 ) -1 ) #prior for Dirichlet
)
m_shell_all <- cstan( file= "models/2_shell_knowledgeinput.stan" , data=dat_shells , chains=3, cores = 3 )
# for (i in 1:nrow(d_shellppl)) {
#   if(d_shellppl$anonymeID[i] %in% rownames(d_shell_k)) {
#     d_shellppl$index_k[i] <- which (rownames(d_shell_k) == d_shellppl$anonymeID[i])
#   } else {
#     d_shellppl$index_k[i] <- NA
#   }
# }
# for (i in 1:nrow(d_trapppl)) {
#   if(d_trapppl$anonymeID[i] %in% rownames(d_trap_k)) {
#     d_trapppl$index_k[i] <- which (rownames(d_trap_k) == d_trapppl$anonymeID[i])
#   } else {
#     d_trapppl$index_k[i] <- NA
#   }
# }

dat_shells <- list(
  #knowledge data
  W = 20,#nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = 50,#ncol(d_knowledge$Y_l),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age[1:20], #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2)[1:20],
  answers = d_knowledge$Y_l[1:20, 1:50], #all answers from freelist
  prior_dirichlet = rep( 2, length (1:56 ) -1 ) #prior for Dirichlet
)
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = ifelse( is.na(d_shellppl$position_k), -999, d_shellppl$position_k) , #make vector of length N where is reported position of individuals for whom we have knowledge in the knowledge list
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  age_int = d_shellppl$age ,
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k), #n of individuals for whom we have knowledge
  Q = ncol(d_shell_k),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_shell_k, #all answers from freelist
  prior_dirichlet = rep( 0.2, length (1:56 ) -1 ) #prior for Dirichlet
)
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),
  M = nrow(d_shells),
  ID_i= d_shells$index_id,
  ID_k = ifelse( is.na(d_shellppl$position_k), -999, d_shellppl$position_k) , #make vector of length N where is reported position of individuals for whom we have knowledge in the knowledge list
  returns = as.numeric(d_shells$returns)/1000,
  age = d_shellppl$age / mean(d_shellppl$age),
  age_int = d_shellppl$age ,
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  #height 
  #grip
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  knowledge_impute = ifelse(is.na(d_shellppl$knowledge), 1, 0), #vector of 0/1 for whether knowledge has to be imputed
  #knowledge data
  W = nrow(d_shell_k_ppl), #n of individuals for whom we have knowledge
  Q = ncol(d_knowledge$Y_l),
  O = 56 , #n of ages for which we want to impute knowledge > let's include said?!
  age_irt = d_shell_k_ppl$age, #[W] vector of ages of ppl in irt
  sex_irt = ifelse(d_shell_k_ppl$sex == "m", 1, 2),
  answers = d_knowledge$Y_l, #all answers from freelist
  prior_dirichlet = rep( 0.2, length (1:56 ) -1 ) #prior for Dirichlet
)

#m_shell_irt <- cstan( file= "models/irt_only_3.stan" , data=dat_shells , chains=3, cores = 3 )
