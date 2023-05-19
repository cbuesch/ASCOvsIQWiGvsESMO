################################################################################
#########            Calculating cutoff Values for ASCO                    #########                 
################################################################################
# Investigating which ESMO / IQWiG_RR / Mod-IQWiG_HR category correspond to 
# which ASCO score

#---------------------- Needed Packages, parameters and functions -------------------------
library(tidyverse)
library(cutpointr)
library(vcd)
library(data.table)
library(doParallel)
library(foreach)

# Number of used cores for analysis
num.cl <- 48

# Loading path
path_ana  <- "" # e.g. ".../Simulations/Ana/"

# Saving path of ASCO cutoff values
path_saving <- "" # e.g. ".../Simulations/"



SvenssonCutoff <- function(df., numericVariable = "ASCO.OS", categoricalVariable){
  cutoff <- rep(NA, length(unique(df.[,categoricalVariable]))-1)
  for (i in 1:(length(unique(df.[,categoricalVariable]))-1)) {
    if(i==1){
      cutoff[i] <-
        uniroot(f = function(x){
          table(df.[,categoricalVariable], exclude = levels(df.[,categoricalVariable])[!levels(df.[,categoricalVariable]) %in% as.character(unique(df.[,categoricalVariable]))])[[i]] - sum(df.[,numericVariable]<=x)
        },
        interval = c(0,100000),tol = .Machine$double.eps^0.35)$root
    }else{
      cutoff[i] <-
        uniroot(f = function(x){
          table(df.[,categoricalVariable], exclude = levels(df.[,categoricalVariable])[!levels(df.[,categoricalVariable]) %in% as.character(unique(df.[,categoricalVariable]))])[[i]] - sum(df.[,numericVariable]<=x & df.[,numericVariable]>cutoff[i-1])
        },
        interval = c(cutoff[i-1],100000), tol = .Machine$double.eps^0.35)$root
    }
  }
  return(cutoff)
}

# --> funktion: mit output des weighted cohen kappas, diese Funktion in optimizer reinstecken
# Warning message: "In sqrt((sum(p * sweep(sweep(W, 1, W %*% colSums(p) * (1 - kw)),  :NaNs wurden erzeugt"
# --> durch schlechte start, wodurch die Mehrfeldertafel (table) sehr schief verteilt ist (alles in einer Kategorie)
#     Es wird aber immer noch ein kappa wert ausgegeben, wordurch optim weiter rechnen kann
WeightedCohensKappa_ESMO_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2], x[3]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catESMO := factor(
    ifelse(ASCO<=cutoffs[1], "1",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "2",
                  ifelse(ASCO>cutoffs[2] & ASCO<=cutoffs[3], "3", "4"))),
    levels = c(1,2,3,4))]
  
  # Compute table
  ASCO_ESMO_table <- df.[, table(ESMO, ASCO.catESMO, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_ESMO_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}
WeightedCohensKappa_IQWiGRR_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catIQWiG_RR := factor(
    ifelse(ASCO<=cutoffs[1], "4",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "5", "6")),
    levels = c(4,5,6))]
  
  # Compute table
  ASCO_IQWiGRR_table <- df.[, table(IQWiG_RR, ASCO.catIQWiG_RR, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_IQWiGRR_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}
WeightedCohensKappa_ModIQWiGHR_ASCO  <- function(x, df.){
  # df. has to be a data.table (so that the data steps are faster)!
  
  # order cutoff values, so that it stays logical
  cutoffs <- sort(c(x[1], x[2]))
  
  # Compute ASCO categories using the cutoff values from "cutoffs"
  df.[, ASCO.catModIQWiG_HR := factor(
    ifelse(ASCO<=cutoffs[1], "4",
           ifelse(ASCO>cutoffs[1] & ASCO<=cutoffs[2], "5", "6")),
    levels = c(4,5,6))]
  
  # Compute table
  ASCO_ModIQWiG_HR_table <- df.[, table(Mod.IQWiG_HR, ASCO.catModIQWiG_HR, exclude=NULL)]
  
  # Calculating weighted Cohen Kappa (with quadratic weights)
  WeightedCohenKappa <- Kappa(ASCO_ModIQWiG_HR_table, weights = "Fleiss-Cohen")$Weighted[[1]]
  
  # Returning 99999 if WeightedCohenKappa is negative, otherwise
  # returning negative WeightedCohenKappa (because the optimizer function "optim"
  # performs minimization)
  return(ifelse(WeightedCohenKappa<0, 99999, -WeightedCohenKappa))
}

#------------------------ Calculating ASCO cutoff values ----------------------------
results_cutpoints <- matrix(
  list(), nrow=7, ncol = 21,
  dimnames = list(
    Scenarios = c("StandardScen", "Scen3aWEIB", "Scen3bGOMP", "AllOfTheAbove", 
                  "Scen4", "Scen2_Overpowered", "Scen2_Underpowered"),
    MethodsForCutOffs = c("ESMO_roc01_1", "ESMO_roc01_12", "ESMO_roc01_123",
                          "ESMO_Svensson1", "ESMO_Svensson12", "ESMO_Svensson123",
                          "ESMO_CohensKappa1","ESMO_CohensKappa12","ESMO_CohensKappa123",

                          "IQWiGRR_roc01_4", "IQWiGRR_roc01_45",
                          "IQWiGRR_Svensson4", "IQWiGRR_Svensson45",
                          "IQWiGRR_CohensKappa4","IQWiGRR_CohensKappa45",
                          
                          "Mod.IQWiGHR_roc01_4", "Mod.IQWiGHR_roc01_45",
                          "Mod.IQWiGHR_Svensson4", "Mod.IQWiGHR_Svensson45",
                          "Mod.IQWiGHR_CohensKappa4","Mod.IQWiGHR_CohensKappa45")
    )
  )

for (i in rownames(results_cutpoints)) {
  ### Loading data
  # Choose correct analyzed data sets
  df.names <- if(i=="StandardScen"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "StandardScen"))
  }else if(i=="Scen3aWEIB"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3aWEIB_Ana"))
  }else if(i=="Scen3bGOMP"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3bGOMP_Ana"))
  }else if(i=="AllOfTheAbove"){
    list.files(path = path_ana, pattern = ".rds") %>%
      discard(function(x) str_detect(x, "Scen2")) %>% 
      discard(function(x) str_detect(x, "Scen4"))
  }else if(i=="Scen4"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen4"))
  }else if(i=="Scen2_Overpowered"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen2")) %>% 
      keep(function(x) as.numeric(str_sub(str_split(x, ".HR.var.")[[1]][2], start=1L, end = 3L))<1)
  }else if(i=="Scen2_Underpowered"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen2")) %>% 
      keep(function(x) as.numeric(str_sub(str_split(x, ".HR.var.")[[1]][2], start=1L, end = 3L))>1)
  }
  
  # Load analyzed data sets, select only significant studies
  # and calculate needed additional variables
  df <- paste0(path_ana, df.names) %>% map_dfr(readRDS) %>% filter(sig==1) %>% 
    mutate(
      ESMO.1vs234 = factor(ifelse(ESMO == 1, "1", "234")),
      ESMO.12vs34 = factor(ifelse(ESMO %in% c(1, 2), "12", "34")),
      ESMO.123vs4 = factor(ifelse(ESMO %in% c(1, 2, 3), "123", "4")),
      ESMO = factor(ESMO),
      
      IQWiG_RR.4vs56 = factor(ifelse(IQWiG_RR == 4, "4", "56")),
      IQWiG_RR.45vs6 = factor(ifelse(IQWiG_RR %in% c(4, 5), "45", "6")),
      IQWiG_RR = factor(IQWiG_RR),
      
      Mod.IQWiG_HR.4vs56 = factor(ifelse(Mod.IQWiG_HR == 4, "4", "56")),
      Mod.IQWiG_HR.45vs6 = factor(ifelse(Mod.IQWiG_HR %in% c(4, 5), "45", "6")),
      Mod.IQWiG_HR = factor(Mod.IQWiG_HR)
    )

  ### Calculating ASCO cutoff values for ESMO
    ## ROC: roc01
    results_cutpoints[i, "ESMO_roc01_1"][[1]] <-
      cutpointr(data = df, x = ASCO, class = ESMO.1vs234,
                pos_class = "234", neg_class = "1", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))

    results_cutpoints[i, "ESMO_roc01_12"][[1]] <-
      cutpointr(data = df,x = ASCO, class = ESMO.12vs34,
                pos_class = "34", neg_class = "12", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))

    results_cutpoints[i, "ESMO_roc01_123"][[1]] <-
      cutpointr(data = df,x = ASCO, class = ESMO.123vs4,
                pos_class = "4", neg_class = "123", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))

    ## Svenssons method
    SvenssonsCutoff_ESMO <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "ESMO")
    results_cutpoints[i, "ESMO_Svensson1"][[1]]   <- SvenssonsCutoff_ESMO[1]
    results_cutpoints[i, "ESMO_Svensson12"][[1]]  <- SvenssonsCutoff_ESMO[2]
    results_cutpoints[i, "ESMO_Svensson123"][[1]] <- SvenssonsCutoff_ESMO[3]

    ## Maximizing cohens kappa approach
    # Local Minima with different start values (--> parallel computing)
    n_start_values <- 40
    set.seed(123456)
    start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 4))
    start_values[,1] <- c(
      c("roc01+Paper", "Svensson+Paper", "roc01", "Svensson", "MeanOfGroups"),
      rep("Random", n_start_values-5))
    start_values[1:5,2:4]  <- matrix(c(
      # roc01 + Paper
      mean(results_cutpoints[i, "ESMO_roc01_1"][[1]]$optimal_cutpoint), 42.250, 44.815,
      # Svensson + Paper
      results_cutpoints[i, "ESMO_Svensson1"][[1]], 42.250, 44.815,
      # roc01
      mean(results_cutpoints[i, "ESMO_roc01_1"][[1]]$optimal_cutpoint),
      mean(results_cutpoints[i, "ESMO_roc01_12"][[1]]$optimal_cutpoint),
      mean(results_cutpoints[i, "ESMO_roc01_123"][[1]]$optimal_cutpoint),
      # Svensson
      results_cutpoints[i, "ESMO_Svensson1"][[1]],
      results_cutpoints[i, "ESMO_Svensson12"][[1]],
      results_cutpoints[i, "ESMO_Svensson123"][[1]],
      # mean of group means
      mean(c(mean(df$ASCO[which(df$ESMO=="1")]), mean(df$ASCO[which(df$ESMO=="2")]))),
      mean(c(mean(df$ASCO[which(df$ESMO=="2")]), mean(df$ASCO[which(df$ESMO=="3")]))),
      mean(c(mean(df$ASCO[which(df$ESMO=="3")]), mean(df$ASCO[which(df$ESMO=="4")])))
    ),
    ncol = 3, byrow = TRUE)

    for (j in 6:n_start_values) {
      start_values[j,2:4] <- sort(runif(n = 3, min = min(df$ASCO), max = max(df$ASCO)))
    }

    # Saving data as "data.table", selecting only needed coloumns and keep df as "data.frame"
    df_help <- setDT(df)[, list(ASCO, ESMO)]
    df <- as.data.frame(df)

    # setting up parallel computing
    cl <- makeCluster(num.cl) # Specify how many cores should be used
    registerDoParallel(cl)

    # starting parallel computing
    MaxKappa_ESMO_ASCO <- 
      foreach(j = iter(start_values, by='row'),
             .packages = c("vcd", "data.table"), .combine = rbind) %dopar% {
               opt_NelderMead <- optim(fn=WeightedCohensKappa_ESMO_ASCO, par = j[2:4],
                                       df. = df_help,
                                       method = "Nelder-Mead",
                                       control = list(
                                         # Relative convergence tolerance:
                                         # The algorithm stops if it is unable to reduce the value
                                         # by a factor of reltol * (abs(val) + reltol) at a step.
                                         # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                         reltol = sqrt(.Machine$double.eps)
                                       ))
               
               return <- data.frame(start_value_name = j[1],
                                    start_value1 = j[2],
                                    start_value2 = j[3],
                                    start_value3 = j[4],
                                    optimal_value1 = sort(opt_NelderMead$par)[1],
                                    optimal_value2 = sort(opt_NelderMead$par)[2],
                                    optimal_value3 = sort(opt_NelderMead$par)[3],
                                    Converge       = opt_NelderMead$convergence,
                                    KappaValue     = opt_NelderMead$value
               )
               colnames(return) <- c("start_values_name", "start_values1", "start_values2", "start_values3",
                                     colnames(return)[5:9])
               
               return(return)
               }
    # stopping parallel computing
    stopCluster(cl)
   
    # Order maximizing Cohens Kappa results / output
    MaxKappa_ESMO_ASCO_orderd <- MaxKappa_ESMO_ASCO[order(MaxKappa_ESMO_ASCO$KappaValue),]

    # get optimal cutoffs and put them into "results_cutpoints"
    results_cutpoints[i, "ESMO_CohensKappa1"][[1]]   <- MaxKappa_ESMO_ASCO_orderd[1,c(5,9)]
    results_cutpoints[i, "ESMO_CohensKappa12"][[1]]  <- MaxKappa_ESMO_ASCO_orderd[1,c(6,9)]
    results_cutpoints[i, "ESMO_CohensKappa123"][[1]] <- MaxKappa_ESMO_ASCO_orderd[1,c(7,9)]

   
  ### Calculating ASCO cutoff values for IQWiGRR
    ## ROC: roc01
    results_cutpoints[i, "IQWiGRR_roc01_4"][[1]] <-
      cutpointr(data = df, x = ASCO, class = IQWiG_RR.4vs56,
                pos_class = "56", neg_class = "4", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))
      
    results_cutpoints[i, "IQWiGRR_roc01_45"][[1]] <-
      cutpointr(data = df,x = ASCO, class = IQWiG_RR.45vs6,
                pos_class = "6", neg_class = "45", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))

    ## Svenssons method
    SvenssonsCutoff_IQWiGRR <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "IQWiG_RR")
    results_cutpoints[i, "IQWiGRR_Svensson4"][[1]]   <- SvenssonsCutoff_IQWiGRR[1]
    results_cutpoints[i, "IQWiGRR_Svensson45"][[1]]  <- SvenssonsCutoff_IQWiGRR[2]
    
    ## Maximizing cohens kappa approach
    # Local Minima with different start values (--> parallel computing)
    n_start_values <- 40
    set.seed(123456)
    start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 3))
    start_values[,1] <- c(
      c("roc01", "Svensson", "MeanOfGroups"),
      rep("Random", n_start_values-3))
    start_values[1:3,2:3]  <- matrix(c(
      # roc01
      mean(results_cutpoints[i, "IQWiGRR_roc01_4"][[1]]$optimal_cutpoint), mean(results_cutpoints[i, "IQWiGRR_roc01_45"][[1]]$optimal_cutpoint),
      # Svensson
      results_cutpoints[i, "IQWiGRR_Svensson4"][[1]], results_cutpoints[i, "IQWiGRR_Svensson45"][[1]],
      # mean of group means
      mean(c(mean(df$ASCO[which(df$IQWiG_RR=="4")]), mean(df$ASCO[which(df$IQWiG_RR=="5")]))),
      mean(c(mean(df$ASCO[which(df$IQWiG_RR=="5")]), mean(df$ASCO[which(df$IQWiG_RR=="6")])))
    ),
    ncol = 2, byrow = TRUE)
    
    for (j in 4:n_start_values) {
      start_values[j,2:3] <- sort(runif(n = 2, min = min(df$ASCO), max = max(df$ASCO)))
    }
    
    # Saving data as "data.table", selecting only needed coloumns and keep df as "data.frame"
    df_help <- setDT(df)[, list(ASCO, IQWiG_RR)]
    df <- as.data.frame(df)
    
    # Setting up parallel computing
    cl <- makeCluster(num.cl) # Specify how many cores should be used
    registerDoParallel(cl)
    
    # Starting parallel computing
    MaxKappa_IQWiGRR_ASCO <- foreach(j = iter(start_values, by='row'),
                                   .packages = c("vcd", "data.table"), .combine = rbind) %dopar% {
                                     opt_NelderMead <- optim(fn=WeightedCohensKappa_IQWiGRR_ASCO, par = j[2:3],
                                                             df. = df_help,
                                                             method = "Nelder-Mead",
                                                             control = list(
                                                               # Relative convergence tolerance:
                                                               # The algorithm stops if it is unable to reduce the value
                                                               # by a factor of reltol * (abs(val) + reltol) at a step.
                                                               # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                                               reltol = sqrt(.Machine$double.eps)#,
                                                               #trace = 6
                                                             ))
                                     return <- data.frame(start_values_name = j[1],
                                                          start_value1 = j[2],
                                                          start_value2 = j[3],
                                                          optimal_value1 = sort(opt_NelderMead$par)[1],
                                                          optimal_value2 = sort(opt_NelderMead$par)[2],
                                                          Converge       = opt_NelderMead$convergence,
                                                          KappaValue     = opt_NelderMead$value)
                                     colnames(return) <- c("start_values_name", "start_values1", "start_values2",
                                                           colnames(return)[4:7])
                                     
                                     return(return)
                                   }
    # Stopping parallel computing
    stopCluster(cl)
    
    # Order maximizing Cohens Kappa results / output
    MaxKappa_IQWiGRR_ASCO_orderd <- MaxKappa_IQWiGRR_ASCO[order(MaxKappa_IQWiGRR_ASCO$KappaValue),]
    
    # Get optimal cutoffs and put them into "results_cutpoints"
    results_cutpoints[i, "IQWiGRR_CohensKappa4"][[1]]  <- MaxKappa_IQWiGRR_ASCO_orderd[1, c(4,7)]
    results_cutpoints[i, "IQWiGRR_CohensKappa45"][[1]] <- MaxKappa_IQWiGRR_ASCO_orderd[1, c(5,7)]
    
    
  ### Calculating ASCO cutoff values for Mod.IQWiGHR
    ## ROC: roc01
    results_cutpoints[i, "Mod.IQWiGHR_roc01_4"][[1]] <-
      cutpointr(data = df, x = ASCO, class = Mod.IQWiG_HR.4vs56,
                pos_class = "56", neg_class = "4", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))
    
    results_cutpoints[i, "Mod.IQWiGHR_roc01_45"][[1]] <-
      cutpointr(data = df,x = ASCO, class = Mod.IQWiG_HR.45vs6,
                pos_class = "6", neg_class = "45", direction = ">=",
                method = minimize_metric, metric = roc01, break_ties = c) %>%
      select(optimal_cutpoint, roc01) %>%
      unnest(cols = c(optimal_cutpoint, roc01))

    ## Svenssons method 
    SvenssonsCutoff_Mod.IQWiGHR <- SvenssonCutoff(df. = as.data.frame(df), numericVariable = "ASCO", categoricalVariable = "Mod.IQWiG_HR")
    results_cutpoints[i, "Mod.IQWiGHR_Svensson4"][[1]]  <- SvenssonsCutoff_Mod.IQWiGHR[1]
    results_cutpoints[i, "Mod.IQWiGHR_Svensson45"][[1]] <- SvenssonsCutoff_Mod.IQWiGHR[2]

    ## Maximizing cohens kappa approach
      # Local Minima with different start values (--> parallel computing)
      n_start_values <- 40
      set.seed(123456)
      start_values <- data.frame(matrix(NA, nrow = n_start_values, ncol = 3))
      start_values[,1] <- c(
        c("roc01", "Svensson", "MeanOfGroups"),
        rep("Random", n_start_values-3))
      start_values[1:3,2:3]  <- matrix(c(
        # roc01
        mean(results_cutpoints[i, "Mod.IQWiGHR_roc01_4"][[1]]$optimal_cutpoint), mean(results_cutpoints[i, "Mod.IQWiGHR_roc01_45"][[1]]$optimal_cutpoint),
        # Svensson
        results_cutpoints[i, "Mod.IQWiGHR_Svensson4"][[1]], results_cutpoints[i, "Mod.IQWiGHR_Svensson45"][[1]],
        # mean of group means
        mean(c(mean(df$ASCO[which(df$Mod.IQWiG_HR=="4")]), mean(df$ASCO[which(df$Mod.IQWiG_HR=="5")]))),
        mean(c(mean(df$ASCO[which(df$Mod.IQWiG_HR=="5")]), mean(df$ASCO[which(df$Mod.IQWiG_HR=="6")])))
        ),
        ncol = 2, byrow = TRUE)
      
      for (j in 4:n_start_values) {
        start_values[j,2:3] <- sort(runif(n = 2, min = min(df$ASCO), max = max(df$ASCO)))
      }

      # Saving data as "data.table", selecting only needed coloumns and keep df itsel as "data.frame"
      df_help <- setDT(df)[, list(ASCO, Mod.IQWiG_HR)]
      df <- as.data.frame(df)
      
      # Setting up parallel computing
      cl <- makeCluster(num.cl) # Specify how many cores should be used
      registerDoParallel(cl)
      
      # Starting parallel computing
      MaxKappa_ModIQWiGHR_ASCO <- foreach(j = iter(start_values, by='row'),
                                        .packages = c("vcd", "data.table"), .combine = rbind) %dopar% {
                                          opt_NelderMead <- optim(fn=WeightedCohensKappa_ModIQWiGHR_ASCO, par = j[2:3],
                                                                  df. = df_help,
                                                                  method = "Nelder-Mead",
                                                                  control = list(
                                                                    # Relative convergence tolerance:
                                                                    # The algorithm stops if it is unable to reduce the value
                                                                    # by a factor of reltol * (abs(val) + reltol) at a step.
                                                                    # Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
                                                                    reltol = sqrt(.Machine$double.eps)#,
                                                                    #trace = 6
                                                                    ))
                                          return <- data.frame(start_values_name = j[1],
                                                               start_value1 = j[2],
                                                               start_value2 = j[3],
                                                               optimal_value1 = sort(opt_NelderMead$par)[1],
                                                               optimal_value2 = sort(opt_NelderMead$par)[2],
                                                               Converge       = opt_NelderMead$convergence,
                                                               KappaValue     = opt_NelderMead$value)
                                          colnames(return) <- c("start_values_name", "start_values1", "start_values2",
                                                                colnames(return)[4:7])
                                          return(return)
                                          }
      # Stopping parallel computing
      stopCluster(cl)

      # Order maximizing Cohens Kappa results / output
      MaxKappa_ModIQWiGHR_ASCO_ordered <- MaxKappa_ModIQWiGHR_ASCO[order(MaxKappa_ModIQWiGHR_ASCO$KappaValue),]

      # Get optimal cutoffs and put them into "results_cutpoints"
      results_cutpoints[i, "Mod.IQWiGHR_CohensKappa4"][[1]]  <- MaxKappa_ModIQWiGHR_ASCO_ordered[1, c(4,7)]
      results_cutpoints[i, "Mod.IQWiGHR_CohensKappa45"][[1]] <- MaxKappa_ModIQWiGHR_ASCO_ordered[1, c(5,7)]
}

### saving the complete results
saveRDS(results_cutpoints,
        file = paste0(path_saving, "Results_ASCOCutoffValues.rds"))


