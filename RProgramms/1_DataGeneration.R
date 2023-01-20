################################################################################
######                           Data generation                          ######
################################################################################

# Setting working directory (where functions and Simulation.seed are saved)
setwd("") 

#---------------------- Needed Packages and Functions --------------------------
library(foreach)
library(doParallel)

source("0_CostumeFunctions_DataGeneration.R")

#-------------------- Same parameter in all scenarios  -------------------------
# Number of used cores for Simulations
num.cl <- 48

# Number of iterations for each sub-scenario
n_sim <- 10000 

# Seed for simulation
load("Simulation.seed.RData")
#seed = sample(x = 1:1000000000, size = it, replace = FALSE)
#save(seed, file = "Simulation.seed.RData")

# Path for saving generated data:
path_data <- "" # e.g. ".../Simulations/Data/"


#-------------------------- Standard Scenario ----------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- 0.6 
accrual.time <- 24 # months --> 2 years
fu.time      <- NA # will be set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, 
                          med.C, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

# Calculating parameter lambda for the independent exponential censoring (lambda.cens)
parameters$lambda.cens <- apply(
  X = parameters, 
  MARGIN = 1, 
  FUN = function(x){
    rootSolve::uniroot.all(
      f = function(aMax, lambda1, p1, lambda2, p2, duration, cens.rate, lambdaC){
        res <- 1/aMax * (
          p1*lambda1/((lambda1+lambdaC)^2) * ( exp((aMax-duration)*(lambda1+lambdaC)) - exp(-duration*(lambda1+lambdaC)) ) +
            p2*lambda2/((lambda2+lambdaC)^2) * ( exp((aMax-duration)*(lambda2+lambdaC)) - exp(-duration*(lambda2+lambdaC)) )
          ) + 
          p1 + p2 - p1*lambda1/(lambda1+lambdaC) - p2*lambda2/(lambda2+lambdaC) - 
          cens.rate
        return(res)
      }, 
      interval  = c(-1,1), 
      aMax      = x["accrual.time"], 
      lambda1   = log(2)/x["med.C"],
      p1        = 0.5,
      lambda2   = (x["HR"]*x["HR.var"])*log(2)/x["med.C"], 
      p2        = 0.5,
      duration  = x["accrual.time"] + x["fu.time"],
      cens.rate = x["cens.rate"])
    }
  )

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters$n.T,decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters, file = paste0(path_data, "StandardScen_parameters.rds"))

#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                       data.scen.i[[i]] <- 
                         DataEXP(
                           n.T            = para$n.T, 
                           n.C            = para$n.C,
                           median.control = para$med.C,
                           accrual.time   = para$accrual.time,
                           fu.time        = para$fu.time,
                           HR             = para$HR, 
                           HR.var         = para$HR.var, 
                           lambda.cens    = para$lambda.cens,
                           index.seed     = seed[i])
                      }
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
                                           "StandardScen_Data", 
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)

#------------------------------ Scenario 2 -------------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- c(0.8, 0.9, 1.1, 1.2)
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- 0.6 
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")

# Calculating parameter lambda for the independent exponential censoring (lambda.cens)
parameters$lambda.cens <- apply(
  X = parameters,
  MARGIN = 1, 
  FUN = function(x){
    rootSolve::uniroot.all(
      f = function(aMax, lambda1, p1, lambda2, p2, duration, cens.rate, lambdaC){
        res <- 1/aMax * (
          p1*lambda1/((lambda1+lambdaC)^2) * ( exp((aMax-duration)*(lambda1+lambdaC)) - exp(-duration*(lambda1+lambdaC)) ) +
            p2*lambda2/((lambda2+lambdaC)^2) * ( exp((aMax-duration)*(lambda2+lambdaC)) - exp(-duration*(lambda2+lambdaC)) )
          ) + 
          p1 + p2 - p1*lambda1/(lambda1+lambdaC) - p2*lambda2/(lambda2+lambdaC) - 
          cens.rate
        return(res)
        }, 
      interval  = c(-1,1), 
      aMax      = x["accrual.time"], 
      lambda1   = log(2)/x["med.C"],
      p1        = 0.5,
      lambda2   = (x["HR"]*x["HR.var"])*log(2)/x["med.C"], 
      p2        = 0.5,
      duration  = x["accrual.time"] + x["fu.time"],
      cens.rate = x["cens.rate"])
    }
  )

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters$n.T, decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen2_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataEXP(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              lambda.cens    = para$lambda.cens,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data,
                                           "Scen2_Data",
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3a: WEIBULL --------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- 0.6
shape.C.T    <- c(0.5, 1.5)   
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later


# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          shape.C.T, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,9],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "shape.C.T", "HR.var", "HR", "n_sim", "n.T", 
                          "n.C")

# Calculating parameter lambda for the independent exponential censoring (lambda.cens)
Integrate_lambda.cens_Weib <- function(aMax, dur, p1, k1, l1, p2, k2, l2, cens.rate, lC){
  non.ad.cens.2 <- 1/aMax * 
    integrate(function(a, dur, p1, k1, l1, p2, k2, l2, lC){
      sapply(a, 
             function(a){
               integrate(function(t, a, dur, p1, k1, l1, p2, k2, l2, lC) {
                 return(
                   p1*k1*l1*(t*l1)^(k1-1) * exp(-(t*l1)^k1) * exp(-lC*t) + 
                     p2*k2*l2*(t*l2)^(k2-1) * exp(-(t*l2)^k2) * exp(-lC*t)
                 )
               }, lower = 0, upper = dur-a,
               a = a, dur = dur, 
               p1 = p1, k1 = k1, l1 = l1, 
               p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
             }
      )
    }, lower = 0, upper = aMax,
    dur = dur, 
    p1 = p1, k1 = k1, l1 = l1, 
    p2 = p2, k2 = k2, l2 = l2, lC = lC)$value
  
  return(p1 + p2 - non.ad.cens.2 - cens.rate)
}

parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = Integrate_lambda.cens_Weib, 
          interval  = c(0,1), extendInt = "yes",
          aMax      = x["accrual.time"], 
          dur       = x["accrual.time"] + x["fu.time"],
          p1        = 0.5,
          k1        = x["shape.C.T"],
          l1        = ((log(2))^(1/x["shape.C.T"])) / x["med.C"],
          p2        = 0.5,
          k2        = x["shape.C.T"],
          l2        = ((x["HR"]*x["HR.var"]*log(2)) / (x["med.C"]^x["shape.C.T"]))^(1/x["shape.C.T"]), 
          cens.rate = x["cens.rate"])$root
  }
)

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,11],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen3aWEIB_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "n_sim" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataWEIB(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              lambda.cens    = para$lambda.cens,
                              shape.C.T      = para$shape.C.T,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
                                           "Scen3aWEIB_Data",
                                           "..beta.",     para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.", para$cens.rate,
                                           ".r.",         para$r,
                                           ".med.C.",     para$med.C,
                                           ".shape.C.T.", para$shape.C.T,
                                           ".HR.var.",    para$HR.var,
                                           ".HR.",        para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------- Scenario 3b: Gompertz -------------------------------
#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- 0.6 
shape.C.T    <- c(-0.2, 0.2)   
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          shape.C.T, HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,9],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", "r", 
                          "med.C", "shape.C.T",  "HR.var", "HR", "n_sim", "n.T", 
                          "n.C")

# Calculating parameter lambda for the independent exponential censoring (lambda.cens)
Integrate_lambda.cens_Gomp <- function(aMax, dur, p1, a1, b1, p2, a2, b2, cens.rate, lC){
  non.ad.cens2 <- 1/aMax * integrate(function(a, dur, p1, a1, b1, p2, a2, b2, lC){
    sapply(a, 
           function(a){
             integrate(function(t, a, dur, p1, a1, b1, p2, a2, b2, lC) {
               return(
                 (
                   p1*b1*exp(a1*t)*exp( -b1/a1*(exp(a1*t)-1) )  +
                     p2*b2*exp(a2*t)*exp( -b2/a2*(exp(a2*t)-1) )
                 ) * exp(-lC*t)
               )
             }, lower = 0, upper = dur-a,
             a = a, dur = dur, 
             p1 = p1, a1 = a1, b1 = b1, 
             p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
           }
    )
  }, lower = 0, upper = aMax,
  dur = dur, 
  p1 = p1, a1 = a1, b1 = b1, 
  p2 = p2, a2 = a2, b2 = b2, lC = lC)$value
  
  return(p1+p2-non.ad.cens2 - cens.rate)
}

parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = Integrate_lambda.cens_Gomp, 
          interval  = c(0,1), extendInt = "yes",
          aMax      = x["accrual.time"], 
          dur       = x["accrual.time"] + x["fu.time"],
          p1        = 0.5,
          a1        = x["shape.C.T"],
          b1        = x["shape.C.T"]*log(2) / (exp(x["med.C"]*x["shape.C.T"])-1),
          p2        = 0.5,
          a2        = x["shape.C.T"],
          b2        = (x["HR"]*x["HR.var"])*x["shape.C.T"]*log(2) / (exp(x["med.C"]*x["shape.C.T"])-1),
          cens.rate = x["cens.rate"])$root
})

# Excluding scenarios where lambda.cens < 0 (administrative censoring already 
# too strong so that a larger censoring rate than wanted is present):
parameters <- parameters[which(parameters$lambda.cens>0),]

# Sorting parameter matrix after the sample size
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,11],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen3bGOMP_parameters.rds"))


#### Data generation
# Starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival", "flexsurv")) %dopar% {
                     # Generating data "it" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataGomp(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              lambda.cens    = para$lambda.cens,
                              shape.C.T      = para$shape.C.T,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data, 
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
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)


#------------------------------ Scenario 4 -------------------------------------
# treat.effect.start is set to 1/3*med.C due to ESMOs disadvantage of using 
# "gain"(med.T-med.C)

#### Parameters
med.C        <- c(6,12,18,24,30) # in months
HR           <- seq(0.3, 0.9, by = 0.02)
HR.var       <- 1
beta         <- c(0.1, 0.2)
alpha        <- 0.05
r            <- 1
cens.rate    <- 0.6
accrual.time <- 24 # in months --> 2 years
fu.time      <- NA # will set to 2*med.C later

# Creating parameter matrix with all sub-scenarios
parameters <- expand.grid(beta, accrual.time, fu.time, cens.rate, r, med.C, 
                          HR.var, HR, n_sim)

# Setting follow up time to 2*med.C later
parameters[,3] <- 2*parameters[,6]

# Calculating sample size for each sub-scenario and adding it to the parameter 
# matrix
parameters$n.T <- NA
parameters$n.C <- NA
for (i in 1:dim(parameters)[1]) {
   n <-
      sample.size.func(
         alpha          = alpha,
         beta           = parameters[i,1],
         n.proportion   = parameters[i,5],
         HR             = parameters[i,8],
         median.control = parameters[i,6],
         fu.time        = parameters[i,3],
         accrual.time   = parameters[i,2],
         cens.rate      = parameters[i,4])
   parameters$n.T[i] = n$n.T
   parameters$n.C[i] = n$n.C
}

# Adding column names to parameter matrix
colnames(parameters) <- c("beta", "accrual.time", "fu.time", "cens.rate", 
                          "r", "med.C", "HR.var", "HR", "n_sim", "n.T", "n.C")


# Calculating parameter lambda for the independent exponential censoring (lambda.cens)
parameters$lambda.cens <- apply(X = parameters, MARGIN = 1, FUN = function(x){
  uniroot(f = function(aMax, dur, startT, p1, l1, p2, l2, cens.rate, lC){
    cens <- 1/aMax * (
      p1*l1/((l1+lC)^2) * ( exp(-(dur-aMax)*(l1+lC)) - exp(-dur*(l1+lC)) )
      + p2*l2/((l2+lC)^2) * exp(-l1*startT) * ( exp(-l2*(dur-aMax-startT) - lC*(dur-aMax)) - exp(-l2*(dur-startT) - lC*dur) )
      ) + p1 - p2*( exp(-l2*startT) - exp(-l1*startT) - 1 ) - p1*l1/(l1+lC) - p2*l2/(l2+lC) * ( exp(-startT*(l1 + lC)) - exp(-startT*(l2+lC)) + 1 )
    return(cens - cens.rate)
    },
    interval  = c(0,1), extendInt = "yes",
    aMax      = x["accrual.time"], 
    dur       = x["accrual.time"] + x["fu.time"],
    startT    = 1/3 * x["med.C"],
    p1        = 0.5,
    l1        = log(2) / x["med.C"],
    p2        = 0.5,
    l2        = (x["HR"]*x["HR.var"])*log(2) / x["med.C"],
    cens.rate = x["cens.rate"])$root
  }
)

# Sorting parameter matrix after the sample size 
# (so that the parallel computing runs faster)
parameters <- parameters[order(parameters[,10],decreasing=TRUE),]

# Saving parameter matrix
saveRDS(parameters,
        file = paste0(path_data, "Scen4_parameters.rds"))

#### Data generation
# starting "timer" for simulation
start.time.1 <- Sys.time()

# Setting up parallel computing
cl <- makeCluster(num.cl)
registerDoParallel(cl)

# Starting parallel computing
results = foreach(para = iter(parameters, by='row'), 
                  .combine = rbind, .packages = c("survival")) %dopar% {
                     # Generating data "it" times for each sub-scenario
                     data.scen.i <- rep(list(NA), para$n_sim)
                     for (i in 1:para$n_sim) {
                        data.scen.i[[i]] <- 
                           DataEXPNonProp(
                              n.T            = para$n.T, 
                              n.C            = para$n.C,
                              median.control = para$med.C,
                              accrual.time   = para$accrual.time,
                              fu.time        = para$fu.time,
                              HR             = para$HR, 
                              HR.var         = para$HR.var, 
                              effect.start.T = 1/3 * para$med.C,
                              lambda.cens    = para$lambda.cens,
                              index.seed     = seed[i])
                     }
                     
                     # Saving data set
                     saveRDS(data.scen.i,
                             file = paste0(path_data,
                                           "Scen4_Data",
                                           "..beta.",        para$beta,
                                           ".accrual.time.", para$accrual.time,
                                           ".fu.time.",      para$fu.time,
                                           ".cens.rate.",    para$cens.rate,
                                           ".r.",            para$r,
                                           ".med.C.",        para$med.C,
                                           ".HR.var.",       para$HR.var,
                                           ".HR.",           para$HR,
                                           ".rds"))
                  }
# Stopping parallel computing 
stopCluster(cl)

# Calculating needed time used for simulation
(time.1 <- Sys.time() - start.time.1)
