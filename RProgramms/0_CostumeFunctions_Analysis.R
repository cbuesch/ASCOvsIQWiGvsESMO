################################################################################
####        Costume functions needed for analysis of generated data         ####
################################################################################

#--------------------- Analysis of generated data sets -------------------------
analysis.assessment.methods <- function(df){
  #### Analysis of generated data sets
  ## Input:
  # df: Generated data set
  
  ## Output: 
  # HR.point:     Estimated hazard ratio point estimate
  # HR.CI.low:    Estimated lower limit of the 95% confidence interval of the 
  #               hazard ratio
  # HR.CI.up:     Estimated upper limit of the 95% confidence interval of the 
  #               hazard ratio
  # median.C:     Median survival time of control group
  # median.T:     Median survival time of treatment group 
  # median.gain:  Difference of median.T - median.C
  # surv.gain:    Survival rate difference between both treatment groups 
  #               (time point depending on observed median.C)
  # ASCO:         ASCO's additional benefit assessment method
  # IQWiG_RR:     IQWiG's additional benefit assessment method
  # Mod.IQWiG_HR: Modified IQWiG's additional benefit assessment method
  # ESMO:         ESMO's additional benefit assessment method
  # sig:          Significant difference (in the correct direction) between 
  #               treatment groups present
  # If data set does not show a significant difference between treatment and 
  # control group using the log-rank test, "sig" is set to 0 (not significant)
  # and all other outputs are set to NA.
  
  ## Information for data set construction
  # Data.TandC: Survival / Censoring time
  # Data.status: Indication of censoring and event/death 
  #              (0: censored, 1: event/death)
  # Data.Z: Treatment allocation / group 
  
  #### Function start
  ### Perform Cox-Regression
  cox.sum <- summary(coxph(Surv(Data.TandC, Data.status != 0) ~ Data.Z, 
                           data = df))
  ### Checking if Treatment is significant in the correct direction
  if(cox.sum$waldtest[[3]]>0.05 | cox.sum$conf.int[3]>1){
    ## Not significant, hence no further analysis is performed
    Ana.out <- data.frame(HR.point = NA, HR.CI.low = NA, HR.CI.up = NA,
                          median.C = NA, median.T = NA, median.gain = NA,
                          surv.gain = NA,
                          ASCO = NA, IQWiG_RR = NA, Mod.IQWiG_HR = NA, ESMO = NA,
                          sig = 0)
  }else{
    ## Significant, hence further analysis is performed
    
    ## Calculation of HR, upper and lower CI
    HR.point  <- cox.sum$conf.int[1]
    HR.CI.low <- cox.sum$conf.int[3]
    HR.CI.up  <- cox.sum$conf.int[4]
    
    ## Calculation of median survival times
    # If the a survival curve of the simulated data did not fall underneath 50%, 
    # a conservative approach was implemented using the last observed maximal 
    # event or censoring time present instead. 
    survtest <- summary(
      survfit(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df)
    )$table
    # Treatment group (median.T)
    if(is.na(survtest[2,7])){
      median.T <- max(df$Data.TandC[df$Data.Z=="Treatment"])
    }else{
      median.T <- survtest[2,7]
    }
    # Control group (median.C)
    if(is.na(survtest[1,7])){
      median.C <- max(df$Data.TandC[df$Data.Z=="Control"])
    }else{
      median.C <- survtest[1,7]
    }
    # median.gain
    median.gain <- abs(median.T - median.C)
    
    ## Calculation of med.C (decides which part of ESMO is used)
    if (median.C <=12){
      med.C = "<=12"
    }else if (median.C > 12 && median.C <= 24) {
      med.C = "12.24"
    } else if (median.C > 24) {
      med.C = ">24"
    }
    
    ## Survival gain
    # Some of the survival curves of the simulated data did not reach 2, 3 or 5
    # year survival and hence the ESMO score could not be calculated. To solve 
    # this issue we extended the survival curve, meaning that either the last
    # censoring time was carried forward or a survival rate of 0% was used.
    fit.C <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, 
                     data = df[which(df$Data.Z=="Control"),])
    fit.T <- survfit(Surv(Data.TandC, Data.status != 0) ~ 1, 
                     data = df[which(df$Data.Z=="Treatment"),])
    
    surv.gain <-
      switch (as.character(med.C),
              "<=12" = {
                surv.C = summary(fit.C, times=2*12, extend=TRUE)$surv 
                surv.T = summary(fit.T, times=2*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              "12.24" = {
                surv.C = summary(fit.C, times=3*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=3*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              },
              ">24" = {
                surv.C = summary(fit.C, times=5*12, extend=TRUE)$surv
                surv.T = summary(fit.T, times=5*12, extend=TRUE)$surv
                if(is.null(surv.C) && is.null(surv.T)){
                  0
                }else if(is.null(surv.C)){
                  abs(0 - surv.T)
                }else if(is.null(surv.T)){
                  abs(surv.C - 0)
                }else{
                  abs(surv.C - surv.T)
                }
              })
    
    ## Additional benefit methods
    # ASCO
      # Clinical Benefit score (CB) 
      ASCO_CB <- (1 - HR.point)*100
      # Bonus Point (BP): Tail of the survival curve 
        # Time point on the survival curve that is 2 times the estimated median 
        # of the control group
        time_point <- median.C*2
        # Proportion of patients alive at time_point
        help <- summary(survfit(Surv(Data.TandC, Data.status != 0) ~ Data.Z, data = df), times = time_point, extend = TRUE)
        p_C  <- help$surv[1]
        p_T  <- help$surv[2]
        # ASCO_BP
        if( (p_C>=0.2) && # assuming > 20% surviving in control group, otherwise no bonus points
            ( (p_T - p_C) / p_C >= 1.5) # 50% or greater improvement in proportion of patients alive with the treatment group
            ){
          ASCO_BP <- 20
          }else{
            ASCO_BP <- 0
          }
      # Net Health Benefit score (NHB)
        ASCO <- ASCO_CB + ASCO_BP
        
    # IQWiG_RR
    IQWiG_RR <- (HR.CI.up<0.85)*6 + 
      (HR.CI.up<0.95 & HR.CI.up>=0.85)*5 +
      (HR.CI.up<1 & HR.CI.up>=0.95)*4
    
    # Mod.IQWiG-HR
    Mod.IQWiG_HR <- (HR.CI.up<0.7908765)*6 + 
      (HR.CI.up<0.9286669 & HR.CI.up>=0.7908765)*5 + 
      (HR.CI.up<1 & HR.CI.up>=0.9286669)*4
    
    # ESMO
    ESMO <-
      switch (as.character(med.C),
              "<=12" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.7 || median.gain<1.5)*1 +
                    (  (HR.CI.low<=0.65 && median.gain<2 && median.gain>=1.5) ||
                         (HR.CI.low<0.7 && HR.CI.low>0.65 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.65 && median.gain<3 && median.gain>=2)*3 +
                    (HR.CI.low<=0.65 && median.gain>=3)*4
                }
              },
              "12.24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<1.5)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<3 && median.gain>=1.5) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=1.5)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<5 && median.gain>=3)*3 +
                    (HR.CI.low<=0.7 && median.gain>=5)*4
                }
              },
              ">24" = {
                if(surv.gain>=0.1){
                  4
                }else{
                  (HR.CI.low>0.75 || median.gain<4)*1 +
                    ( (HR.CI.low<=0.7 && median.gain<6 && median.gain>4) ||
                        (HR.CI.low<0.75 && HR.CI.low>0.7 && median.gain>=4)  )*2 +
                    (HR.CI.low<=0.7 && median.gain<9 && median.gain>=6)*3 +
                    (HR.CI.low<=0.7 && median.gain>=9)*4
                }
              })
    
    ## Analyzed data
    Ana.out <- data.frame(HR.point, HR.CI.low, HR.CI.up,
                          median.C, median.T, median.gain,
                          surv.gain,
                          ASCO, IQWiG_RR, Mod.IQWiG_HR, ESMO,
                          sig = 1)
  }
  ### Return analyzed data 
  return(Ana.out)
}