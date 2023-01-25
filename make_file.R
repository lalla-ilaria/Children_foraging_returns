
rm(list = ls())

source("init_project.R")

source("project_support.R")

tic()
if (dir.exists("raw_data")) source("0_clean_data.R")
toc()

tic()
source("1_simulation.R") # ~25 mins
toc()

tic()
source("2_analysis.R")  # ~3 hours
toc()

tic()
source("3_plots.R")
toc()

tic()
source("4_supp_info_plots.R") 
toc()
