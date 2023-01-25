
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

source("raw_data/corrections.R")

d_trips        <- read.csv("raw_data/trips_sheet.csv")
d_activity     <- read.csv("raw_data/activity_sheet.csv")
d_returns      <- read.csv("raw_data/returns_sheet.csv")
d_participants <- read.csv("raw_data/participants_sheet.csv")

d_trips        <- correct_trips(d_trips)
d_activity     <- correct_activity(d_activity)
d_returns      <- correct_returns(d_returns)
d_participants <- correct_participants(d_participants)

d_anthropometrics <- read.csv("raw_data/anthropometrics.csv")

census <- read.csv("raw_data/census_BK_2020.csv")

d_knowledge <- list.load("raw_data/processed_data.RData")

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
dat_shells <- make_list_data_age(foraging_type = "shells")
dat_traps <- make_list_data_age(foraging_type = "traps")

post_s <- extract.samples(m_shell_age)
post_t <- extract.samples(m_trap_age)

#generated data

#proportion max foraging by age 10
p_max_for_10 <- matrix(NA, 
                       nrow = 3, ncol = 2, 
                       dimnames = list(c("median", "5%PI", "94%PI"), 
                                       c("shells","traps")))
p_max_for_10[1,1] <- median((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ median(post_s$gamma))
p_max_for_10[2:3,1] <- PI(((1-exp(-post_s$beta * 10/mean_age_shells  )) ^ post_s$gamma))

p_max_for_10[1,2] <- (1-exp(-median(post_t$beta) * 10/mean_age_traps  )) ^ median(post_t$gamma)
p_max_for_10[2:3,2] <- PI((1-exp(-post_t$beta * 10/mean_age_traps  )) ^ PI(post_t$gamma))
#WHY MEDIAN IS NOT BETWEEN PI????




######################################
# TIDE LEVELS
######################################
dat_shells <- make_list_data_all(foraging_type = "shells")

dat_tides <- dat_shells [ c("M", "ID_i", "tide", "age", "sex")]
dat_tides$age <- d_tides$age[d_tides$ID_i]
dat_tides$sex <- d_tides$sex[d_tides$ID_i]
dat_tides$age <- dat_tides$age * mean_age_shells

dat_tides <- as.data.frame(dat_tides)
dat_tides <- dat_tides [-which(dat_tides$age >=15 & dat_tides$sex == "m"),]
dat_tides %>% filter(age >=20 ) %>% summarize( Mean = mean(tide))
dat_tides %>% filter(age <=19 ) %>% summarize( Mean = mean(tide))

dat_tides %>% filter(tide >0 ) %>% summarize( Mean = mean(age))
dat_tides %>% filter(tide <=0 ) %>% summarize( Mean = mean(age))

###############
#shells calculation
dat_shells <- make_list_data_all(foraging_type = "shells")
post_s <- extract.samples(m_shells_all)

phi <-  apply(post_s$iota,1,mean )  +
  post_s$gamma * log(1-exp(- post_s$beta * 1.23  )) +
  post_s$eta_h* mean(log(dat_shells$height), na.rm = TRUE) +
  post_s$theta_g* mean(log(dat_shells$grip), na.rm = TRUE) +
  post_s$zeta_k* min(apply(post_s$knowledge, 2, mean))
psi <-  post_s$xi * mean(log(180/mean(real_data$shells$lenght_min))) +
  post_s$tau* 0
kg_shells <- exp (log(post_s$alpha) + phi_min + psi +
                (post_s$sigma^2 /2))
mean(kg_shells)
PI(kg_shells)

##########
#SAVE DATA
##########
generated_data <- list(p_max_for_10 = p_max_for_10)

list.save(generated_data, "3_outcomes/generated_data.Rdata")
