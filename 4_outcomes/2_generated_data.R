library(rethinking)
library(rlist)
library(tidyverse)
library(ggridges)
real_data <- list.load("2_data_preparation/processed_data.RData")
seq_trait <- seq(0,3,0.001)


anonyme <- function( df_column) {
  anonyme_id <- vector( length = length(df_column) )
  for(i in 1:length(df_column)){
    anonyme_id[i] <- paste( which(letters == tolower(substring(df_column[i], 1, 1))), parse_number(df_column[i]), sep = "")
  }
  return(anonyme_id)
}


#########################
#DESCRIPTIVE STATISTICS 
#for M&M - non reproducible
#########################

source("../data/corrections.R")

d_trips        <- read.csv("../data/trips_sheet.csv")
d_activity     <- read.csv("../data/activity_sheet.csv")
d_returns      <- read.csv("../data/returns_sheet.csv")
d_participants <- read.csv("../data/participants_sheet.csv")

d_trips        <- correct_trips(d_trips)
d_activity     <- correct_activity(d_activity)
d_returns      <- correct_returns(d_returns)
d_participants <- correct_participants(d_participants)

d_anthropometrics <- read.csv("../data/anthropometrics.csv")

census <- read.csv("../data/census_BK_2020.csv")

d_knowledge <- list.load("../data/processed_data.RData")

#n and percentage of hunting trips
sum(str_detect(d_trips$trip_type, "fyuka"))/nrow(d_trips)

#n households
hhs <- unique(census$household_2019)
hhs <- hhs[- which (hhs == "BK144, NY201" |
                    hhs == "BK109, BK106" |
                    hhs == "BK134, BK161" |
                    hhs == "BK121, BK142" |
                    hhs == "inferred" |
                    hhs == "out_of_village" |
                    hhs == "BK144, BK172, BK176" |
                    hhs == "BK158, BK159" |
                    hhs == "dead" |
                    hhs == "REMOVED" |
                    hhs == "" |
                    is.na(hhs)
)]
length(hhs)

#people in village
nrow(census) - sum(census$household_2019 %in% c("inferred", "out_of_village", "dead", "REMOVED", ""))

#sample sizes
d_knowledge$N #people in the knowledge sample
sum(!is.na(d_anthropometrics$height)) #people for whom we have height
sum(!is.na(d_anthropometrics$weight)) #people for whom we have weight
sum(!is.na(d_anthropometrics$grip)) #people for whom we have grip strength

#######
#descriptive shellfish foraging
shells <- d_returns %>% 
  filter( trip_type %in% c("shell_collect", "crab_collect") ) %>% 
  select(trip, who,weight_g )
#remove trips where weight has not been measured
shells <- shells %>%
  filter( !is.na(weight_g))
length(unique(shells$trip)) #number of shellfish foraging trips 
nrow(shells) #number of foraging trips/person 
length(unique(shells$who))

#descriptive sample
shell_ppl <- data.frame(ID = unique(shells$who),
                        age = rep(NA, length( unique(shells$who))),
                        sex = rep(NA, length( unique(shells$who)))
)
for(i in 1:nrow(shell_ppl)){
  shell_ppl$age[i] <- census$age[ which (census$ID == shell_ppl$ID[i])]
}
for(i in 1:nrow(shell_ppl)){
  shell_ppl$sex[i] <- census$sex[ which (census$ID == shell_ppl$ID[i])]
}
#proportion females
sum(shell_ppl$sex == "f")/nrow(shell_ppl)
#age range
min(shell_ppl$age, na.rm = T)
max(shell_ppl$age, na.rm = T)
mean(shell_ppl$age, na.rm = T)

#n trips per age group
for(i in 1:nrow(shells)){
  shells$age[i] <- shell_ppl$age[ which (shell_ppl$ID == shells$who[i])]
}

sum(shells$age <=19, na.rm = T)
sum(shells$age >=20, na.rm = T)


#######
trap_trips <- d_trips$trip[ which (str_detect(d_trips$trip_type, "fyuka"))]
#n people who participated to trap trips
length(d_participants$ID[which(d_participants$trip %in% trap_trips)])#total trips/participant
length(unique(d_participants$ID[which(d_participants$trip %in% trap_trips)]))#n unique participants
#n trips
length(trap_trips)

#other measures in wrangle data

#Traps
people_traps <- unique(d_participants$ID[which(d_participants$trip %in% trap_trips)])
trap_ppl <- data.frame(ID = people_traps,
                       age = rep(NA, length( people_traps )),
                       sex = rep(NA, length( people_traps ))
)
trap_ppl$anonymeID <- anonyme(trap_ppl$ID)

for(i in 1:nrow(trap_ppl)){
  trap_ppl$age[i] <- census$age[ which (census$ID == trap_ppl$ID[i])]
}

for(i in 1:nrow(trap_ppl)){
  trap_ppl$sex[i] <- census$sex[ which (census$ID == trap_ppl$ID[i])]
}

trap_ppl$success <- ifelse(trap_ppl$anonymeID %in% real_data$traps$anonymeID, 1, 0)

#age range - all participants
min(trap_ppl$age, na.rm = T)
max(trap_ppl$age, na.rm = T)
mean(trap_ppl$age, na.rm = T)

#age range -only who set traps
min(trap_ppl$age[which(trap_ppl$success == 1)], na.rm = T)
max(trap_ppl$age[which(trap_ppl$success == 1)], na.rm = T)
mean(trap_ppl$age[which(trap_ppl$success == 1)], na.rm = T)

#n traps captured something
length(unique(real_data$individual_traps$trap_ID))#tot n traps
length(unique(real_data$individual_traps$trap_ID[which(real_data$individual_traps$success == 1)]))#traps that captured
#total n success
sum(real_data$traps$success)

#group sizes
#function to get mode of distribution
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

n_members_shells <- shells %>% group_by(trip) %>% count()
mean(n_members_shells$n)
median(n_members_shells$n)
getmode(n_members_shells$n)

n_members_traps <- d_participants %>% filter(trip %in% trap_trips) %>% group_by(trip) %>% count()
mean(n_members_traps$n)
median(n_members_traps$n)
getmode(n_members_traps$n)

#########################
#AGE ONLY
#########################
#prepare data age only
dc_shellppl <- real_data$shell_ppl[complete.cases(real_data$shell_ppl$age),]
dc_shell_k <- real_data$shell_k[which(rownames(real_data$shell_k) %in% dc_shellppl$anonymeID),]
dc_shells <- real_data$shells[which(real_data$shells$anonymeID %in% dc_shellppl$anonymeID),]

dc_shellppl$index_id <- as.integer(as.factor(dc_shellppl$anonymeID))
dc_shellppl <- dc_shellppl[order(dc_shellppl$index_id),]
dc_shells$index_id <- as.integer(as.factor(dc_shells$anonymeID))
dc_shell_k <- dc_shell_k[ order(as.factor(row.names(dc_shell_k))), ]

dc_trapppl <- real_data$trap_ppl[complete.cases(real_data$trap_ppl$age),]
dc_traps <- real_data$traps[which(real_data$traps$anonymeID %in% dc_trapppl$anonymeID),]
dc_trap_k <- real_data$trap_k[which(rownames(real_data$trap_k) %in% dc_trapppl$anonymeID),]

dc_trapppl$index_id <- as.integer(as.factor(dc_trapppl$anonymeID))
dc_trapppl <- dc_trapppl[order(dc_trapppl$index_id),]
dc_traps$index_id <- as.integer(as.factor(dc_traps$anonymeID))
dc_trap_k <- dc_trap_k[ order(as.factor(row.names(dc_trap_k))), ]

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

#generated data

#proportion max foraging by age 10
p_max_for_10 <- matrix(NA, 
                       nrow = 3, ncol = 2, 
                       dimnames = list(c("median", "5%PI", "94%PI"), 
                                       c("shells","traps")))
p_max_for_10[1,1] <- median((1-exp(-post_s$beta * 10/mean(d_shellppl$age)  )) ^ median(post_s$gamma))
p_max_for_10[2:3,1] <- PI(((1-exp(-post_s$beta * 10/mean(d_shellppl$age)  )) ^ post_s$gamma))

p_max_for_10[1,2] <- (1-exp(-median(post_t$beta) * 10/mean(d_trapppl$age)  )) ^ median(post_t$gamma)
p_max_for_10[2:3,2] <- PI((1-exp(-post_t$beta * 10/mean(d_trapppl$age)  )) ^ PI(post_t$gamma))
#WHY MEDIAN IS NOT BETWEEN PI????




######################################
# TIDE LEVELS
######################################
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells

for(i in 1:nrow(d_shells)){
  d_shells$age[i] <- d_shellppl$age[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}

for(i in 1:nrow(d_shells)){
  d_shells$sex[i] <- d_shellppl$sex[which (d_shellppl$anonymeID == d_shells$anonymeID[i])]
}
d_tides <- d_shells [-which(d_shells$age >=15 & d_shells$sex == "m"),]
d_tides %>% filter(age >=20 ) %>% summarize( Mean = mean(tide_avg_depth))
d_tides %>% filter(age <=19 ) %>% summarize( Mean = mean(tide_avg_depth))

d_tides %>% filter(tide_avg_depth >0 ) %>% summarize( Mean = mean(age))
d_tides %>% filter(tide_avg_depth <=0 ) %>% summarize( Mean = mean(age))

###############
#shells calculation

#SHELLS
d_shellppl <- real_data$shell_ppl
d_shells <- real_data$shells
d_shell_k <- real_data$shell_k

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

#SHELLS
dat_shells <- list(
  #foraging data
  N = nrow(d_shellppl),                       #n individuals in total sample
  M = nrow(d_shells),                         #n trip/person
  ID_i= d_shells$index_id,                    #index of person of trip 
  returns = as.numeric(d_shells$returns)/1000,#amount of shells in kg
  age = (d_shellppl$age / mean(d_shellppl$age)),
  sex = ifelse(d_shellppl$sex == "m", 1, 2), #make vector of sexes 1 = male 2 = female
  duration = d_shells$lenght_min/mean(d_shells$lenght_min),
  tide = d_shells$tide_avg_depth,
  #height data
  has_height = ifelse(is.na(d_shellppl$height), 0, 1),# #vector of 0/1 for whether height has to be imputed
  height = d_shellppl$height/mean(d_shellppl$height, na.rm = TRUE),
  min_height = 50/mean(d_shellppl$height, na.rm = TRUE),#average height of newborn as intercept in height model
  #grip data
  has_grip = ifelse(is.na(d_shellppl$grip), 0, 1),# #vector of 0/1 for whether grip has to be imputed
  grip = d_shellppl$grip/mean(d_shellppl$grip, na.rm = TRUE),
  #knowledge data
  has_knowledge = ifelse(is.na(d_shellppl$knowledge), 0, 1),# #vector of 0/1 for whether knowledge has to be imputed
  Q = ncol(d_shell_k),                        #n items in freelist
  answers = d_shell_k                         #all answers from freelist
)

phi <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(180/mean(d_shells$lenght_min))) +
  post_s$tau* 0
kg_shells <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
mean(kg_shells)
PI(kg_shells)

##########
#SAVE DATA
##########
generated_data <- list(p_max_for_10 = p_max_for_10)

list.save(generated_data, "4_outcomes/generated_data.Rdata")
