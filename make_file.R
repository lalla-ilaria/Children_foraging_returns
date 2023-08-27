##############################
######## MAKEFILE ############
##############################

#####################
#prepare environment#
#####################

#set working directory to "Children_foraging_returns/"

#init project 
#remove folders and files generated during compilation of the project
if (dir.exists("plots")) unlink("plots", recursive = TRUE)
if (dir.exists("4_outcomes/model_fit")) unlink("4_outcomes/model_fit", recursive = TRUE)

model_binaries <- list.files("models", full.names = TRUE)
model_binaries <- model_binaries[-grep("\\.stan$", model_binaries)]
file.remove(model_binaries)

#create folders necessary in the directory
dir.create("plots")
dir.create("plots/validate_model")
dir.create("4_outcomes/model_fit")

#load packages

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

usePackage("dagitty")
usePackage("rlist")
usePackage("tidyverse")
usePackage("ggnewscale")
usePackage("truncnorm")
usePackage("cowplot")

library(rethinking)
library(dagitty)
library(rlist)
library(tidyverse)
library(ggnewscale)#for using two scales of color with ggplot
library(truncnorm)#for truncated normal
library(cowplot)#for pasting plots together

#define colors and other plotting things
real_data <- list.load("2_data_preparation/processed_data.RData")
tide_data <- read.csv("2_data_preparation/tide_data.csv")
boycol  <- rgb(114/255, 181/255, 40/255) #"navyblue"
girlcol <- rgb(208/255, 29/255, 157/255) #"red3"
trapcol <- "#52b502" #"#5ca81e"
shellcol <-  "#ea5914"
othercol <- "grey30"##1482ac"
seq_trait <- seq(0,3,0.02)
age_plot <- 40
      
options(bitmapType='cairo')

#add line to get clean slate file.delete and tell which files

##############
##simulation##
##############
#compile models and create plots - simulated data
#input: simulation code
#output: plots and values for model validation
#time: 32.52 min
start_time <- Sys.time()
source("1_simulation/1_simulation.R")
source("1_simulation/2_analysis.R")
end_time <- Sys.time()
time_1data <- end_time - start_time
time_1data




###############
####analyses###
###############
#input:
#output:
#time 
#loads function to prepare data
source("2_data_preparation/2_prep_data_functions.R")
#compile models and create plots - real data
start_time <- Sys.time()
source("3_analysis/1_analysis.R")
source("4_outcomes/1_plots.R")
source("4_outcomes/2_generated_data.R")
source("4_outcomes/3_supp_info_plots.R")
end_time <- Sys.time()
time_2data <- end_time - start_time
time_2data

