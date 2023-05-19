################################################################################
#########            Analysis of the generated data                    #########                 
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(tidyverse)
library(foreach)
library(doParallel)

source("0_CostumeFunctions_Analysis.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for analysis
num.cl <- 48

# Loading / Saving paths
path_data <- "" # e.g. ".../Simulations/Data/"
path_ana  <- "" # e.g. ".../Simulations/Ana/"


#------------------------- Standard Scenario 1 ---------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "StandardScen_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_StandardScen <- foreach(para = iter(parameters, by='row'), 
                        .packages = "survival", .combine = rbind) %dopar% {
   # Load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(
     paste0(
       path_data,
       "StandardScen_Data",
       "..beta.",        para$beta,
       ".accrual.time.", para$accrual.time,
       ".fu.time.",      para$fu.time,
       ".cens.rate.",    para$cens.rate,
       ".r.",            para$r,
       ".med.C.",        para$med.C,
       ".HR.var.",       para$HR.var,
       ".HR.",           para$HR,
       ".rds"
       )
     )
   
   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, 
                                                  ncol = dim(ana.data.scen.1)[2],
                                                  nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)
   
   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Adding parameters to Scenario:
   rownames(para) <- c()
   results.i <- cbind(para, ana.data.scen.i)
   
   # Saving of analysis results
   saveRDS(results.i,
           file = paste0(
             path_ana, 
             "StandardScen_Ana",
             "..beta.",        para$beta,
             ".accrual.time.", para$accrual.time,
             ".fu.time.",      para$fu.time,
             ".cens.rate.",    para$cens.rate,
             ".r.",            para$r,
             ".med.C.",        para$med.C,
             ".HR.var.",       para$HR.var,
             ".HR.",           para$HR,
             ".rds"
             )
           )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 2 -------------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen2_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_Scen2 <- foreach(para = iter(parameters, by='row'), 
                  .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_times" times replicates
   df.loaded <- readRDS(
     paste0(
       path_data,
       "Scen2_Data",
       "..beta.",        para$beta,
       ".accrual.time.", para$accrual.time,
       ".fu.time.",      para$fu.time,
       ".cens.rate.",    para$cens.rate,
       ".r.",            para$r,
       ".med.C.",        para$med.C,
       ".HR.var.",       para$HR.var,
       ".HR.",           para$HR,
       ".rds"
       )
     )

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, 
                                                  ncol = dim(ana.data.scen.1)[2], 
                                                  nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Adding parameters to Scenario:
   rownames(para) <- c()
   results.i <- cbind(para, ana.data.scen.i)
   
   # Saving of analysis results
   saveRDS(results.i,
           file = paste0(
             path_ana,
             "Scen2_Ana",
             "..beta.",        para$beta,
             ".accrual.time.", para$accrual.time,
             ".fu.time.",      para$fu.time,
             ".cens.rate.",    para$cens.rate,
             ".r.",            para$r,
             ".med.C.",        para$med.C,
             ".HR.var.",       para$HR.var,
             ".HR.",           para$HR,
             ".rds"
             )
           )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3a: WEIBULL --------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen3aWEIB_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results_Scen3a <- foreach(para = iter(parameters, by='row'),
                        .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(
     paste0(
       path_data,
       "Scen3aWEIB_Data",
       "..beta.",        para$beta,
       ".accrual.time.", para$accrual.time,
       ".fu.time.",      para$fu.time,
       ".cens.rate.",    para$cens.rate,
       ".r.",            para$r,
       ".med.C.",        para$med.C,
       ".shape.C.T.",    para$shape.C.T,
       ".HR.var.",       para$HR.var,
       ".HR.",           para$HR,
       ".rds"
       )
     )

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Adding parameters to Scenario:
   rownames(para) <- c()
   results.i <- cbind(para, ana.data.scen.i)
   
   # Saving of analysis results
   saveRDS(results.i,
           file = paste0(
             path_ana, 
             "Scen3aWEIB_Ana",
             "..beta.",        para$beta,
             ".accrual.time.", para$accrual.time,
             ".fu.time.",      para$fu.time,
             ".cens.rate.",    para$cens.rate,
             ".r.",            para$r,
             ".med.C.",        para$med.C,
             ".shape.C.T.",    para$shape.C.T,
             ".HR.var.",       para$HR.var,
             ".HR.",           para$HR,
             ".rds"
             )
           )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3b: GOMPERTZ -------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen3bGOMP_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results_Scen3b <- foreach(para = iter(parameters, by='row'), 
                         .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(
     paste0(
       path_data,
       "Scen3bGOMP_Data",
       "..beta.",        para$beta,
       ".accrual.time.", para$accrual.time,
       ".fu.time.",      para$fu.time,
       ".cens.rate.",    para$cens.rate,
       ".r.",            para$r,
       ".med.C.",        para$med.C,
       ".shape.C.T.",    para$shape.C.T,
       ".HR.var.",       para$HR.var,
       ".HR.",           para$HR,
       ".rds"
       )
     )
   
   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Adding parameters to Scenario:
   rownames(para) <- c()
   results.i <- cbind(para, ana.data.scen.i)
   
   # Saving of analysis results
   saveRDS(results.i,
           file = paste0(
             path_ana, 
             "Scen3bGOMP_Ana",
             "..beta.",        para$beta,
             ".accrual.time.", para$accrual.time,
             ".fu.time.",      para$fu.time,
             ".cens.rate.",    para$cens.rate,
             ".r.",            para$r,
             ".med.C.",        para$med.C,
             ".shape.C.T.",    para$shape.C.T,
             ".HR.var.",       para$HR.var,
             ".HR.",           para$HR,
             ".rds"
             )
           )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 4 -------------------------------------
## Loading Parameters
parameters  <- readRDS(paste0(path_data, "Scen4_parameters.rds"))

## Analysis
# Starting "timer" for analysis
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results_Scen4 <- foreach(para = iter(parameters, by='row'), 
                  .packages = "survival", .combine = rbind) %dopar% {
   # load data of sub-scenarios with "n_sim" times replicates
   df.loaded <- readRDS(
     paste0(
       path_data,
       "Scen4_Data",
       "..beta.",        para$beta,
       ".accrual.time.", para$accrual.time,
       ".fu.time.",      para$fu.time,
       ".cens.rate.",    para$cens.rate,
       ".r.",            para$r,
       ".med.C.",        para$med.C,
       ".HR.var.",       para$HR.var,
       ".HR.",           para$HR,
       ".rds"
       )
     )

   # Analysis of each of the "n_sim" times simulated trial
   ana.data.scen.1           <- analysis.assessment.methods(df = df.loaded[[1]])
   ana.data.scen.i           <- data.frame(matrix(NA, ncol = dim(ana.data.scen.1)[2], nrow = para$n_sim))
   ana.data.scen.i[1,]       <- ana.data.scen.1
   colnames(ana.data.scen.i) <- names(ana.data.scen.1)

   for (i in 2:para$n_sim) {
      ana.data.scen.i[i,] <- analysis.assessment.methods(df = df.loaded[[i]])
   }
   
   # Adding parameters to Scenario:
   rownames(para) <- c()
   results.i <- cbind(para, ana.data.scen.i)
   
   # Saving of analysis results
   saveRDS(results.i,
           file = paste0(
             path_ana, 
             "Scen4_Ana",
             "..beta.",        para$beta,
             ".accrual.time.", para$accrual.time,
             ".fu.time.",      para$fu.time,
             ".cens.rate.",    para$cens.rate,
             ".r.",            para$r,
             ".med.C.",        para$med.C,
             ".HR.var.",       para$HR.var,
             ".HR.",           para$HR,
             ".rds"
             )
           )
}
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for analysis
(time.1 <- Sys.time() - start.time.1)