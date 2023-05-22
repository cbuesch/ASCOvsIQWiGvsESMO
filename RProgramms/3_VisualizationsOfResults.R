################################################################################
######                  Visualization of analyzed data                    ######
################################################################################

# Setting working directory 
path_ana  <- "" # e.g. ".../Simulations/Ana/"
path_data <- "" # e.g. ".../Simulations/Data/"
path_saving_cutoffs <- "" # Saving path of ASCO cutoff values, e.g. ".../Simulations/"


#----------------- Needed Packages and further setting ups ---------------------
library(tidyverse)
library(patchwork)  # for convenient plotting
library(ggpubr)     # for convenient plotting
library(flexsurv)   # for calculating hazards of Gompertz distribution
library(pcaPP)      # to speed up kendall tau estimations
library(foreach)    # for parallel computing
library(doParallel) # for parallel computing

#------------- Figure 2: Spearman correlation Standard Scenario ----------------
## Loading results from "2_DataAnalysis_ASCO_CutoffValues.R"
results_cutpoints <- readRDS(paste0(path_saving_cutoffs, "Results_ASCOCutoffValues.rds")) 

## New ASCO categorical variables from cutoff points
# ESMO
ESMO_CohensKappa_Exp <- sort(round(
  c(results_cutpoints["StandardScen", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["StandardScen", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["StandardScen", "ESMO_CohensKappa123"][[1]][[1]]),
  0))
ESMO_CohensKappa_Weib <- sort(round(
  c(results_cutpoints["Scen3aWEIB", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["Scen3aWEIB", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["Scen3aWEIB", "ESMO_CohensKappa123"][[1]][[1]]),
  0))
ESMO_CohensKappa_Gomp <- sort(round(
  c(results_cutpoints["Scen3bGOMP", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["Scen3bGOMP", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["Scen3bGOMP", "ESMO_CohensKappa123"][[1]][[1]]),
  0))
ESMO_CohensKappa_NonPropHazards <- sort(round(
  c(results_cutpoints["Scen4", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["Scen4", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["Scen4", "ESMO_CohensKappa123"][[1]][[1]]),
  0))
ESMO_CohensKappa_BiasedHROverpowered <- sort(round(
  c(results_cutpoints["Scen2_Overpowered", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["Scen2_Overpowered", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["Scen2_Overpowered", "ESMO_CohensKappa123"][[1]][[1]]),
  0))
ESMO_CohensKappa_BiasedHRUnderpowered <- sort(round(
  c(results_cutpoints["Scen2_Underpowered", "ESMO_CohensKappa1"][[1]][[1]],
    results_cutpoints["Scen2_Underpowered", "ESMO_CohensKappa12"][[1]][[1]],
    results_cutpoints["Scen2_Underpowered", "ESMO_CohensKappa123"][[1]][[1]]),
  0))

# IQWiGRR
IQWiGRR_CohensKappa_Exp <- sort(round(
  c(results_cutpoints["StandardScen", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["StandardScen", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))
IQWiGRR_CohensKappa_Weib <- sort(round(
  c(results_cutpoints["Scen3aWEIB", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen3aWEIB", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))
IQWiGRR_CohensKappa_Gomp <- sort(round(
  c(results_cutpoints["Scen3bGOMP", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen3bGOMP", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))
IQWiGRR_CohensKappa_NonPropHazards <- sort(round(
  c(results_cutpoints["Scen4", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen4", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))
IQWiGRR_CohensKappa_BiasedHROverpowered <- sort(round(
  c(results_cutpoints["Scen2_Overpowered", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen2_Overpowered", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))
IQWiGRR_CohensKappa_BiasedHRUnderpowered <- sort(round(
  c(results_cutpoints["Scen2_Underpowered", "IQWiGRR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen2_Underpowered", "IQWiGRR_CohensKappa45"][[1]][[1]]),
  0))

# Mod-IQWiGHR
Mod.IQWiGHR_CohensKappa_Exp <- sort(round(
  c(results_cutpoints["StandardScen", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["StandardScen", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))
Mod.IQWiGHR_CohensKappa_Weib <- sort(round(
  c(results_cutpoints["Scen3aWEIB", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen3aWEIB", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))
Mod.IQWiGHR_CohensKappa_Gomp <- sort(round(
  c(results_cutpoints["Scen3bGOMP", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen3bGOMP", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))
Mod.IQWiGHR_CohensKappa_NonPropHazards <- sort(round(
  c(results_cutpoints["Scen4", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen4", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))
Mod.IQWiGHR_CohensKappa_BiasedHROverpowered <- sort(round(
  c(results_cutpoints["Scen2_Overpowered", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen2_Overpowered", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))
Mod.IQWiGHR_CohensKappa_BiasedHRUnderpowered <- sort(round(
  c(results_cutpoints["Scen2_Underpowered", "Mod.IQWiGHR_CohensKappa4"][[1]][[1]],
    results_cutpoints["Scen2_Underpowered", "Mod.IQWiGHR_CohensKappa45"][[1]][[1]]),
  0))

## Start loading data for each scenario and generate plots
plots_ESMO_ASCO <- rep(list(NA), 6)
plots_IQWiGRR_ASCO <- rep(list(NA), 6)
plots_ModIQWiGHR_ASCO <- rep(list(NA), 6)

for (i in c("StandardScen", "Scen3aWEIB", "Scen3bGOMP", "Scen4", "Scen2_Overpowered", "Scen2_Underpowered")) {
  ### Title for subplots/panels
  title_i <- case_when(
    i=="StandardScen"       ~ "Standard Scenario: Exponential",
    i=="Scen3aWEIB"         ~ "Scenario 3a: Weibull",
    i=="Scen3bGOMP"         ~ "Scenario 3b: Gompertz",
    i=="Scen4"              ~ "Scenario 4: Non-proportional hazards",
    i=="Scen2_Overpowered"  ~ "Scenario 2: Overpowered",
    i=="Scen2_Underpowered" ~ "Scenario 2: Underpowered",
    TRUE ~ NA_character_
  )

  line_thickness <- 0.1

  ### Loading data and plotting subplots/panels
  if(i=="StandardScen"){
    df <- as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("StandardScen")) == "StandardScen") %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[1]] <- df %>%
      reframe(perc = 100*prop.table(table(ESMO))) %>% 
      mutate(ESMO = factor(c(1,2,3,4))) %>%
      left_join(df %>%
                  group_by(ESMO) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(ESMO = factor(c(1,2,3,4))),
                by = c("ESMO")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, ESMO = "all")
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>% ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_Exp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Exp[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Exp[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_Exp,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "ESMO", y = "", title = title_i) +
      theme_bw()
    
    plots_IQWiGRR_ASCO[[1]] <- df %>%
      reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      mutate(IQWiG_RR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(IQWiG_RR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(IQWiG_RR = factor(c(4,5,6))),
                by = c("IQWiG_RR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, IQWiG_RR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Exp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Exp[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_Exp,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = expression("IQWiG"[RR]), y = "", title = title_i) +
      theme_bw()
    
    plots_ModIQWiGHR_ASCO[[1]] <- df %>%
      reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>%
      mutate(Mod.IQWiG_HR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(Mod.IQWiG_HR = factor(c(4,5,6))),
                by = c("Mod.IQWiG_HR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, Mod.IQWiG_HR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Exp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Exp[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_Exp,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = expression("Mod-IQWiG"[HR]), y = "", title = title_i) +
      theme_bw() 
    
  }else if(i=="Scen3aWEIB"){
    df <- as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("Scen3aWEIB")) == "Scen3aWEIB") %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[2]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(ESMO))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        ESMO = factor(rep(c(1,2,3,4), times = 3)),
        hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
      ) %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, ESMO) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    ESMO = factor(rep(c(1,2,3,4), times = 3)),
                    hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("ESMO", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            ESMO = "all",
            hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_Weib[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Weib[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Weib[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_Weib,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "ESMO", y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
    plots_IQWiGRR_ASCO[[2]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        IQWiG_RR = factor(rep(c(4,5,6), times = 3)),
        hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
      )  %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, IQWiG_RR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    IQWiG_RR = factor(rep(c(4,5,6), times = 3)),
                    hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("IQWiG_RR", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            IQWiG_RR = "all",
            hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Weib[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Weib[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_Weib,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = expression("IQWiG"[RR]), y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
    plots_ModIQWiGHR_ASCO[[2]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        Mod.IQWiG_HR = factor(rep(c(4,5,6), times = 3)),
        hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
      )  %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    Mod.IQWiG_HR = factor(rep(c(4,5,6), times = 3)),
                    hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("Mod.IQWiG_HR", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            Mod.IQWiG_HR = "all",
            hazard = factor(hazard, level = c("0.5", "1.5", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Weib[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Weib[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_Weib,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = expression("Mod-IQWiG"[HR]), y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
  }else if(i=="Scen3bGOMP"){
    df <- as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("Scen3bGOMP")) == "Scen3bGOMP") %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[3]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(ESMO))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        ESMO = factor(rep(c(1,2,3,4), times = 3)),
        hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
      )  %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, ESMO) %>%
                  summarise(
                    ASCO0 = min(ASCO),
                    ASCO25 = quantile(ASCO, 0.25),
                    ASCO50 = median(ASCO),
                    ASCO75 = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    ESMO = factor(rep(c(1,2,3,4), times = 3)),
                    hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("ESMO", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            ESMO = "all",
            hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_Gomp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Gomp[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_Gomp[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_Gomp,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
    plots_IQWiGRR_ASCO[[3]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        IQWiG_RR = factor(rep(c(4,5,6), times = 3)),
        hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
      ) %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, IQWiG_RR) %>%
                  summarise(
                    ASCO0 = min(ASCO),
                    ASCO25 = quantile(ASCO, 0.25),
                    ASCO50 = median(ASCO),
                    ASCO75 = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    IQWiG_RR = factor(rep(c(4,5,6), times = 3)),
                    hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("IQWiG_RR", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            IQWiG_RR = "all",
            hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Gomp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_Gomp[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_Gomp,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
    plots_ModIQWiGHR_ASCO[[3]] <- df %>%
      rbind(df %>% mutate(shape.C.T="all")) %>%
      group_by(shape.C.T) %>% reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>% 
      rename(hazard = shape.C.T) %>%
      mutate(
        Mod.IQWiG_HR = factor(rep(c(4,5,6), times = 3)),
        hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
      )  %>%
      left_join(df %>%
                  rbind(df %>% mutate(shape.C.T="all")) %>%
                  group_by(shape.C.T, Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0 = min(ASCO),
                    ASCO25 = quantile(ASCO, 0.25),
                    ASCO50 = median(ASCO),
                    ASCO75 = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>%
                  rename(hazard = shape.C.T) %>%
                  mutate(
                    Mod.IQWiG_HR = factor(rep(c(4,5,6), times = 3)),
                    hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
                  ),
                by = c("Mod.IQWiG_HR", "hazard")) %>%
      rbind(
        df %>% rbind(df %>% mutate(shape.C.T="all")) %>% group_by(shape.C.T) %>%
          summarise(
            ASCO0   = min(ASCO),
            ASCO25  = quantile(ASCO, 0.25),
            ASCO50  = median(ASCO),
            ASCO75  = quantile(ASCO, 0.75),
            ASCO100 = max(ASCO)
          ) %>%
          rename(hazard = shape.C.T) %>%
          mutate(
            perc = 100,
            Mod.IQWiG_HR = "all",
            hazard = factor(hazard, level = c("-0.2", "0.2", "all"), labels = c("decreasing", "increasing", "combined"))
          )
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label, fill = hazard)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Gomp[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_Gomp[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_Gomp,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i, fill = "Scenario 3 hazards") +
      theme_bw() 
    
  }else if(i=="Scen4"){
    df <- as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("Scen4")) == "Scen4") %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[4]] <- df %>%
      reframe(perc = 100*prop.table(table(ESMO))) %>% 
      mutate(ESMO = factor(c(1,2,3,4))) %>%
      left_join(df %>%
                  group_by(ESMO) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(ESMO = factor(c(1,2,3,4))),
                by = c("ESMO")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, ESMO = "all")
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_NonPropHazards[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_NonPropHazards[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_NonPropHazards[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_NonPropHazards,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i) +
      theme_bw() 
    
    plots_IQWiGRR_ASCO[[4]] <- df %>%
      reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      mutate(IQWiG_RR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(IQWiG_RR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(IQWiG_RR = factor(c(4,5,6))),
                by = c("IQWiG_RR"))  %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, IQWiG_RR = "all")
      )%>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_NonPropHazards[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_NonPropHazards[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_NonPropHazards,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i) +
      theme_bw() 
    
    plots_ModIQWiGHR_ASCO[[4]] <- df %>%
      reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>% 
      mutate(Mod.IQWiG_HR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(Mod.IQWiG_HR = factor(c(4,5,6))) ,
                by = c("Mod.IQWiG_HR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, Mod.IQWiG_HR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_NonPropHazards[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_NonPropHazards[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_NonPropHazards,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "ASCO", title = title_i) +
      theme_bw() 
    
  }else if(i=="Scen2_Overpowered"){
    df <- as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("Scen2")) == "Scen2") %>% 
      filter(as.numeric(substr(value, 
                               nchar(sub("\\HR.var.*", "", value)) + 8, 
                               nchar(sub("\\HR.var.*", "", value)) + 10)) < 1
      ) %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[5]] <- df %>%
      reframe(perc = 100*prop.table(table(ESMO))) %>% 
      mutate(ESMO = factor(c(1,2,3,4))) %>%
      left_join(df %>%
                  group_by(ESMO) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(ESMO = factor(c(1,2,3,4))),
                by = c("ESMO")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, ESMO = "all")
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHROverpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHROverpowered[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHROverpowered[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_BiasedHROverpowered,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "", y = "", title = title_i) +
      theme_bw() 
    
    plots_IQWiGRR_ASCO[[5]] <- df %>%
      reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      mutate(IQWiG_RR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(IQWiG_RR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(IQWiG_RR = factor(c(4,5,6))),
                by = c("IQWiG_RR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, IQWiG_RR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_BiasedHROverpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_BiasedHROverpowered[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_BiasedHROverpowered,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "", title = title_i)+
      theme_bw()
    
    plots_ModIQWiGHR_ASCO[[5]] <- df %>%
      reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>% 
      mutate(Mod.IQWiG_HR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(Mod.IQWiG_HR = factor(c(4,5,6))),
                by = c("Mod.IQWiG_HR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, Mod.IQWiG_HR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
        
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_BiasedHROverpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_BiasedHROverpowered[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_BiasedHROverpowered,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "", title = title_i) +
      theme_bw() 
    
  }else if(i=="Scen2_Underpowered"){
    as_tibble(list.files(path = path_ana, pattern = ".rds")) %>% 
      filter(substr(value, 1, nchar("Scen2")) == "Scen2") %>% 
      filter(as.numeric(substr(value, 
                               nchar(sub("\\HR.var.*", "", value)) + 8, 
                               nchar(sub("\\HR.var.*", "", value)) + 10)) > 1
      ) %>% 
      sapply(., function(x) paste0(path_ana, x)) %>% 
      map(readRDS) %>% 
      bind_rows() %>% 
      filter(sig == 1)
    
    plots_ESMO_ASCO[[6]] <- df %>%
      reframe(perc = 100*prop.table(table(ESMO))) %>% 
      mutate(ESMO = factor(c(1,2,3,4))) %>%
      left_join(df %>%
                  group_by(ESMO) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(ESMO = factor(c(1,2,3,4))) ,
                by = c("ESMO")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, ESMO = "all")
      ) %>%
      mutate(
        Label = factor(paste0(ESMO, " (", round(perc, 2), "%)"),
                       levels = paste0(ESMO, " (", round(perc, 2), "%)"),
                       labels = paste0(ESMO, " (", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHRUnderpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHRUnderpowered[2]), lwd = line_thickness) +
      geom_hline(aes(yintercept = ESMO_CohensKappa_BiasedHRUnderpowered[3]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=ESMO_CohensKappa_BiasedHRUnderpowered,
                                           name="",
                                           guide = guide_axis(n.dodge=2))) +
      coord_flip() +
      labs(x = "", y = "", title = title_i) +
      theme_bw() 
    
    plots_IQWiGRR_ASCO[[6]] <- df %>%
      reframe(perc = 100*prop.table(table(IQWiG_RR))) %>% 
      mutate(IQWiG_RR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(IQWiG_RR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(IQWiG_RR = factor(c(4,5,6))),
                by = c("IQWiG_RR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, IQWiG_RR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(IQWiG_RR)== "4" ~ "Minor",
          as.character(IQWiG_RR)== "5" ~ "Considerable",
          as.character(IQWiG_RR)== "6" ~ "Major",
          as.character(IQWiG_RR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_BiasedHRUnderpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = IQWiGRR_CohensKappa_BiasedHRUnderpowered[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=IQWiGRR_CohensKappa_BiasedHRUnderpowered,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "", title = title_i) +
      theme_bw()
    
    plots_ModIQWiGHR_ASCO[[6]] <- df %>%
      reframe(perc = 100*prop.table(table(Mod.IQWiG_HR))) %>% 
      mutate(Mod.IQWiG_HR = factor(c(4,5,6))) %>%
      left_join(df %>%
                  group_by(Mod.IQWiG_HR) %>%
                  summarise(
                    ASCO0   = min(ASCO),
                    ASCO25  = quantile(ASCO, 0.25),
                    ASCO50  = median(ASCO),
                    ASCO75  = quantile(ASCO, 0.75),
                    ASCO100 = max(ASCO)
                  ) %>% ungroup() %>% mutate(Mod.IQWiG_HR = factor(c(4,5,6))),
                by = c("Mod.IQWiG_HR")) %>%
      rbind(df %>%
              summarise(
                ASCO0   = min(ASCO),
                ASCO25  = quantile(ASCO, 0.25),
                ASCO50  = median(ASCO),
                ASCO75  = quantile(ASCO, 0.75),
                ASCO100 = max(ASCO)
              ) %>%
              mutate(perc = 100, Mod.IQWiG_HR = "all")
      ) %>%
      mutate(
        help_label = case_when(
          as.character(Mod.IQWiG_HR)== "4" ~ "Minor",
          as.character(Mod.IQWiG_HR)== "5" ~ "Considerable",
          as.character(Mod.IQWiG_HR)== "6" ~ "Major",
          as.character(Mod.IQWiG_HR)== "all" ~ "All",
          TRUE ~ NA_character_
        ),
        Label = factor(paste0(help_label, "\n(", round(perc, 2), "%)"),
                       levels = paste0(help_label, "\n(", round(perc, 2), "%)"),
                       labels = paste0(help_label, "\n(", round(perc, 2), "%)"))
      ) %>%
      ggplot(aes(x = Label)) +
      geom_boxplot(
        aes(ymin = ASCO0, lower = ASCO25, middle = ASCO50, upper = ASCO75, ymax = ASCO100),
        stat = "identity", lwd = line_thickness) +
      # Horizontal lines for ASCO cutoff values
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_BiasedHRUnderpowered[1]), lwd = line_thickness) +
      geom_hline(aes(yintercept = Mod.IQWiGHR_CohensKappa_BiasedHRUnderpowered[2]), lwd = line_thickness) +
      scale_y_continuous(breaks = seq(0,120,by=20), limits = c(0,120),
                         sec.axis=sec_axis(~., breaks=Mod.IQWiGHR_CohensKappa_BiasedHRUnderpowered,
                                           name="",
                                           guide = guide_axis()#guide_axis(n.dodge=2)
                         )
      ) +
      coord_flip() +
      labs(x = "", y = "", title = title_i) +
      theme_bw() 
  }  
}  

## Combination of plots
pFinal_ESMO_ASCO <- plots_ESMO_ASCO[[1]] + plots_ESMO_ASCO[[5]] + plots_ESMO_ASCO[[6]] +
  plots_ESMO_ASCO[[2]] + plots_ESMO_ASCO[[3]] + plots_ESMO_ASCO[[4]] +
  plot_layout(nrow = 2, ncol = 3, guides = "collect") &
  theme(
    plot.margin = margin(0,1,0,0),
    legend.margin=margin(t=-10)
  ) & theme(legend.position = 'bottom')
pFinal_IQWiGRR_ASCO <- plots_IQWiGRR_ASCO[[1]] + plots_IQWiGRR_ASCO[[5]] + plots_IQWiGRR_ASCO[[6]] +
  plots_IQWiGRR_ASCO[[2]] + plots_IQWiGRR_ASCO[[3]] + plots_IQWiGRR_ASCO[[4]] +
  plot_layout(nrow = 2, ncol = 3, guides = "collect") &
  theme(
    plot.margin = margin(0,1,0,0),
    legend.margin=margin(t=-10)
  ) & theme(legend.position = 'bottom')
pFinal_ModIQWiGHR_ASCO <- plots_ModIQWiGHR_ASCO[[1]] + plots_ModIQWiGHR_ASCO[[5]] + plots_ModIQWiGHR_ASCO[[6]] +
  plots_ModIQWiGHR_ASCO[[2]] + plots_ModIQWiGHR_ASCO[[3]] + plots_ModIQWiGHR_ASCO[[4]] +
  plot_layout(nrow = 2, ncol = 3, guides = "collect") &
  theme(
    plot.margin = margin(0,1,0,0),
    legend.margin=margin(t=-10)
  ) & theme(legend.position = 'bottom')


#------------- Figure 5: Spearman correlation Standard Scenario ----------------
# Loading results Standard Scenario and calculating Spearman correlation
Spearman_results <- readRDS(paste0(path_data, "StandardScen_parameters.rds")) %>% 
  mutate(
    Cor.ASCO.IQWiG_RR = NA, Cor.ASCO.ModIQWiG_HR = NA, Cor.ASCO.ESMO = NA
  )

# Setting up parallel computing
num.cl <- 3
cl <- makeCluster(num.cl) 
registerDoParallel(cl)

# Starting parallel computing
Spearman_results <- foreach(para = iter(Spearman_results, by='row'),
                            .combine = rbind, .packages = c("tidyverse")) %dopar% {
                              df.ana <- readRDS(paste0(path_ana, "StandardScen_Ana",
                                                       "..beta.",        para$beta,
                                                       ".accrual.time.", para$accrual.time,
                                                       ".fu.time.",      para$fu.time,
                                                       ".cens.rate.",    para$cens.rate,
                                                       ".r.",            para$r,
                                                       ".med.C.",        para$med.C,
                                                       ".HR.var.",       para$HR.var,
                                                       ".HR.",           para$HR,
                                                       ".rds")) %>%
                                # Using only significant results
                                filter(sig == 1)
                              
                              # Calculating Spearman correlations 
                              # (setting correlation to "n.a." if only the same score was given)
                              para$Cor.ASCO.IQWiG_RR <- ifelse(
                                length(levels(factor(df.ana$IQWiG_RR))) == 1, 
                                "n.a.", 
                                cor(df.ana$ASCO, df.ana$IQWiG_RR, method = "spearman")
                              )
                              para$Cor.ASCO.ModIQWiG_HR <- ifelse(
                                length(levels(factor(df.ana$Mod.IQWiG_HR))) == 1, 
                                "n.a.", 
                                cor(df.ana$ASCO, df.ana$Mod.IQWiG_HR, method = "spearman")
                              )
                              para$Cor.ASCO.ESMO <- ifelse(
                                length(levels(factor(df.ana$ESMO))) == 1, 
                                "n.a.", 
                                cor(df.ana$ASCO, df.ana$ESMO, method = "spearman")
                              )
                              return(para)
                            }
# Stopping parallel computing 
stopCluster(cl)

# Plotting Spearman correlation
Figure5 <- Spearman_results %>% filter(beta == 0.1) %>%
  arrange(beta, cens.rate, med.C, HR) %>%
  mutate(Cor.ASCO.ESMO         = as.numeric(Cor.ASCO.ESMO),
         Cor.ASCO.IQWiG_RR     = as.numeric(Cor.ASCO.IQWiG_RR),
         Cor.ASCO.ModIQWiG_HR  = as.numeric(Cor.ASCO.ModIQWiG_HR),
         
         med.C.text = paste0("Median Control: ", med.C),
         pow.text = paste0(100*(1-as.numeric(as.character(beta))),"% Power"),
         text_together = factor(paste0(pow.text, ", ", med.C.text),
                                levels = c("90% Power, Median Control: 6",  "90% Power, Median Control: 12", "90% Power, Median Control: 18",
                                           "90% Power, Median Control: 24", "90% Power, Median Control: 30"),
                                labels = c("90% Power, Median Control: 6",  "90% Power, Median Control: 12", "90% Power, Median Control: 18",
                                           "90% Power, Median Control: 24", "90% Power, Median Control: 30")
         ),
         HR = as.numeric(as.character(HR))
  ) %>%
  pivot_longer(cols = c(Cor.ASCO.ESMO, Cor.ASCO.IQWiG_RR, Cor.ASCO.ModIQWiG_HR), 
               values_to = "SpearmanCorrelation", 
               names_to = "Methods") %>%
  ggplot(aes(x = HR, y = SpearmanCorrelation, group = Methods, linetype = Methods, shape = Methods)) +
  geom_point(size = 0.2) + geom_line(size = 0.2) + ylim(c(0,1)) +
  labs(title = "", # "Pairwise Spearman correlation between ASCO's and ESMO's / IQWiG’s additional benefit assessment methods",
       x = "designHR", y = "Spearman correlation") +
  scale_linetype_manual(
    name = "Spearman correlation of ",
    values = c(2,4,6),
    labels = c("ASCO vs. ESMO", expression("ASCO vs. IQWiG"[RR]), expression("ASCO vs. Mod-IQWiG"[HR]))
  ) +
  scale_shape_manual(
    name = "Spearman correlation of ",
    values = c(15,16,17),
    labels = c("ASCO vs. ESMO", expression("ASCO vs. IQWiG"[RR]), expression("ASCO vs. Mod-IQWiG"[HR]))
  ) +
  # Setting tick intervals on x-axis
  scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.06), minor_breaks = seq(0.3, 0.9, by = 0.02)) +
  # Setting tick intervals on y-axis
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75), minor_breaks = c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = "bottom"
  ) +
  facet_wrap(~text_together, nrow=1, scales = "fixed")


#--------- Figure 6: Spearman/Kendall tau correlation All Scenarios ------------
results_spearman_kendall_ASCO <-
  as.data.frame(expand.grid(
    Scenario = c("StandardScen",
                 "Scen3aWEIB", "Scen3aWEIBIncreasing", "Scen3aWEIBDecreasing",
                 "Scen3bGOMP", "Scen3bGOMPIncreasing", "Scen3bGOMPDecreasing", 
                 "Scen4", "Scen2_Overpowered", "Scen2_Underpowered"),
    med.C   = c(6, 12, 18, 24, 30, "all"),
    Methods = c("ESMO", "IQWiG_RR", "Mod.IQWiG_HR"))) %>%
  mutate(
    Spearman = NA,
    Kendall  = NA
  )

# Combining results of sub-scenarios and calculating Spearman/Kendall tau correlation 
for(i in levels(results_spearman_kendall_ASCO$Scenario)){
  # Choose correct analyzed data sets
  df.names <- if(i=="StandardScen"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "StandardScen"))
  }else if(i=="Scen3aWEIB"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3aWEIB_Ana"))
  }else if(i=="Scen3aWEIBIncreasing"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3aWEIB_Ana")) %>% 
      keep(function(x) as.numeric(str_sub(str_split(x, ".shape.C.T.")[[1]][2], start=1L, end = 3L))>1)
  }else if(i=="Scen3aWEIBDecreasing"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3aWEIB_Ana")) %>% 
      keep(function(x) as.numeric(str_sub(str_split(x, ".shape.C.T.")[[1]][2], start=1L, end = 3L))<1)
  }else if(i=="Scen3bGOMP"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3bGOMP_Ana"))
  }else if(i=="Scen3bGOMPIncreasing"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3bGOMP_Ana")) %>% 
      keep(function(x) as.numeric(str_split(str_split(x, ".shape.C.T.")[[1]][2], ".HR.var")[[1]][1])>0)
  }else if(i=="Scen3bGOMPDecreasing"){
    list.files(path = path_ana, pattern = ".rds") %>%
      keep(function(x) str_detect(x, "Scen3bGOMP_Ana")) %>% 
      keep(function(x) as.numeric(str_split(str_split(x, ".shape.C.T.")[[1]][2], ".HR.var")[[1]][1])<0)
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
  
  # Load analyzed data sets and select only significant studies
  df <- df.names %>% 
    sapply(., function(x) paste0(path_ana, x)) %>% 
    map(readRDS) %>% 
    bind_rows() %>% 
    filter(sig == 1)
  
  # Calculating Spearman and Kendall correlation 
  for(j in levels(results_spearman_kendall_ASCO$med.C)){
    if(j == "all"){
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "ESMO"), 
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO, df$ESMO, method="spearman"),
          cor.fk(df$ASCO, df$ESMO))
      
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "IQWiG_RR"), 
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO, df$IQWiG_RR, method="spearman"),
          cor.fk(df$ASCO,  df$IQWiG_RR))
      
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "Mod.IQWiG_HR"),
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO,df$Mod.IQWiG_HR, method="spearman"),
          cor.fk(df$ASCO, df$Mod.IQWiG_HR))
    }else{
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "ESMO"),
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO[as.character(df$med.C)==j], df$ESMO[as.character(df$med.C)==j], method="spearman"),
          cor.fk(df$ASCO[as.character(df$med.C)==j], df$ESMO[as.character(df$med.C)==j]))
      
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "IQWiG_RR"),
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO[as.character(df$med.C)==j], df$IQWiG_RR[as.character(df$med.C)==j], method="spearman"),
          cor.fk(df$ASCO[as.character(df$med.C)==j], df$IQWiG_RR[as.character(df$med.C)==j]))
      
      results_spearman_kendall_ASCO[which(as.character(results_spearman_kendall_ASCO$Scenario) == i &
                                            as.character(results_spearman_kendall_ASCO$med.C) == j &
                                            as.character(results_spearman_kendall_ASCO$Methods) == "Mod.IQWiG_HR"),
                                    c("Spearman", "Kendall")] <-
        c(cor(df$ASCO[as.character(df$med.C)==j], df$Mod.IQWiG_HR[as.character(df$med.C)==j], method="spearman"),
          cor.fk(df$ASCO[as.character(df$med.C)==j], df$Mod.IQWiG_HR[as.character(df$med.C)==j]))
      }
  }
}


## Plotting Spearman/Kendall tau correlation All Scenarios using heatmaps
# Spearman correlation
p_ASCO_Spearman <-
  as_tibble(results_spearman_kendall_ASCO) %>% filter(med.C == "all") %>%
  mutate(
    Scenario = factor(Scenario,
                      levels = c("Scen2_Underpowered", "Scen2_Overpowered",
                                 "Scen4",
                                 "Scen3bGOMP", "Scen3bGOMPIncreasing", "Scen3bGOMPDecreasing", 
                                 "Scen3aWEIB", "Scen3aWEIBIncreasing", "Scen3aWEIBDecreasing",
                                 "StandardScen"),
                      labels = c("Scenario 2:\nunderpowered", "Scenario 2:\noverpowered",
                                 "Scenario 4:\nnon-proportional\nhazards", 
                                 "Scenario 3b:\nGompertz", "Scenario 3b:\nincreasing hazards\n(shape > 0)", "Scenario 3b:\ndecreasing hazards\n(shape < 0)",
                                 "Scenario 3a:\nWeibull", "Scenario 3a:\nincreasing hazards\n(shape > 1)", "Scenario 3a:\ndecreasing hazards\n(shape < 1)",
                                 "Standard Scenario:\nExponential"))
  ) %>%
  ggplot(aes(x = Methods, y = Scenario, fill = Spearman))+
  geom_tile() +
  geom_text(aes(label = round(Spearman, 2)), show.legend = FALSE, col = "white") +
  scale_x_discrete("ASCO vs.", labels = c("ESMO", expression(IQWiG[RR]), expression("Mod-IQWiG"[HR]))) +
  labs(fill = "Spearman\ncorrelation") +
  theme_minimal()
p_ASCO_Spearman

# Spearman and Kendall tau correlation 
p_ASCO_SpearmanKendall <-
  as_tibble(results_spearman_kendall_ASCO) %>% filter(med.C == "all") %>%
  pivot_longer(cols = c(Spearman, Kendall), names_to = "AgreementMeasure", values_to = "Values") %>%
  mutate(
    Scenario = factor(Scenario,
                      levels = c("Scen2_Underpowered", "Scen2_Overpowered",
                                 "Scen4",
                                 "Scen3bGOMP", "Scen3bGOMPIncreasing", "Scen3bGOMPDecreasing", 
                                 "Scen3aWEIB", "Scen3aWEIBIncreasing", "Scen3aWEIBDecreasing",
                                 "StandardScen"),
                      labels = c("Scenario 2:\nunderpowered", "Scenario 2:\noverpowered",
                                 "Scenario 4:\nnon-proportional\nhazards",
                                 "Scenario 3b:\nGompertz", "Scenario 3b:\nincreasing hazards\n(shape > 0)", "Scenario 3b:\ndecreasing hazards\n(shape < 0)",
                                 "Scenario 3a:\nWeibull", "Scenario 3a:\nincreasing hazards\n(shape > 1)", "Scenario 3a:\ndecreasing hazards\n(shape < 1)",
                                 "Standard Scenario:\nExponential"))
  ) %>%
  ggplot(aes(x = Methods, y = Scenario, fill = Values))+
  geom_tile() +
  geom_text(aes(label = round(Values, 2)), show.legend = FALSE, col = "white") + 
  scale_x_discrete("ASCO vs.", labels = c("ESMO", expression(IQWiG[RR]), expression("Mod-IQWiG"[HR]))) +
  labs(fill = "Spearman\ncorrelation") +
  theme_minimal() +
  facet_wrap(~AgreementMeasure)
p_ASCO_SpearmanKendall

# Spearman and Kendall tau correlation splitted for med.C
p_ASCO_SpearmanKendall_splitted_medC <-
  as_tibble(results_spearman_kendall_ASCO) %>%
  pivot_longer(cols = c(Spearman, Kendall), names_to = "AgreementMeasure", values_to = "Values") %>%
  mutate(
    Scenario = factor(Scenario,
                      levels = c("Scen2_Underpowered", "Scen2_Overpowered",
                                 "Scen4",
                                 "Scen3bGOMP", "Scen3bGOMPIncreasing", "Scen3bGOMPDecreasing", 
                                 "Scen3aWEIB", "Scen3aWEIBIncreasing", "Scen3aWEIBDecreasing",
                                 "StandardScen"),
                      labels = c("Scenario 2:\nunderpowered", "Scenario 2:\noverpowered",
                                 "Scenario 4:\nnon-proportional\nhazards",
                                 "Scenario 3b:\nGompertz", "Scenario 3b:\nincreasing hazards\n(shape > 0)", "Scenario 3b:\ndecreasing hazards\n(shape < 0)",
                                 "Scenario 3a:\nWeibull", "Scenario 3a:\nincreasing hazards\n(shape > 1)", "Scenario 3a:\ndecreasing hazards\n(shape < 1)",
                                 "Standard Scenario:\nExponential")
                      )
    ) %>%
  ggplot(aes(x = Methods, y = Scenario, fill = Values))+
  geom_tile() +
  geom_text(aes(label = round(Values, 2)), show.legend = FALSE, col = "white") +
  scale_x_discrete("ASCO vs.", labels = c("ESMO", expression(IQWiG[RR]), expression("Mod-IQWiG"[HR]))) +
  labs(fill = "Spearman\ncorrelation") +
  theme_minimal() +
  facet_wrap(med.C~AgreementMeasure)
p_ASCO_SpearmanKendall_splitted_medC


#----------- APPENDIX: Hazard functions for exponential, Weibull ---------------
#---------------- and Gompertz distribution for Appendix -----------------------
# Functions needed for calculation of hazard function
Hazard_Weib_C <- function(shape.C.T.WEIB, time_points, median.control){
   return(
      shape.C.T.WEIB *
         median.control/((log(2))^(1/shape.C.T.WEIB)) *
         (time_points*median.control/((log(2))^(1/shape.C.T.WEIB)))^(shape.C.T.WEIB-1)
   )
}
Hazard_Weib_T <- function(shape.C.T.WEIB, time_points, median.control, HR, HR.var){
   return(
      shape.C.T.WEIB *
         (median.control^shape.C.T.WEIB / (HR*HR.var*log(2)))^(1/shape.C.T.WEIB) *
         (time_points*(median.control^shape.C.T.WEIB / (HR*HR.var*log(2)))^(1/shape.C.T.WEIB))^(shape.C.T.WEIB-1)
   )
}
Hazard_Gomp_C <- function(shape.C.T.GOMP, time_points, median.control){
   return(
      hgompertz(
         x=time_points,
         shape = shape.C.T.GOMP,
         rate = shape.C.T.GOMP*log(2) / (exp(median.control*shape.C.T.GOMP)-1)
      )
   )
}
Hazard_Gomp_T <- function(shape.C.T.GOMP, time_points, median.control, HR, HR.var){
   return(
      hgompertz(
         x=time_points,
         shape = shape.C.T.GOMP,
         rate = (HR*HR.var)*shape.C.T.GOMP*log(2) / (exp(median.control*shape.C.T.GOMP)-1)
      )
   )
}

# Assumed parameters
HR <- 0.9
HR.var <- 1
median.control <- 6
time_points <- seq(1,
                   40,
                   by = 0.2)

# Data frame for plot
df_plot <- data.frame(
   t            = rep(time_points, times=2),
   grp          = c(rep("Control", length(time_points)),
                    rep("Treatment", length(time_points))),
   exp          = c(rep(log(2)/median.control, length(time_points)),
                    rep((HR*HR.var)*log(2)/median.control, length(time_points))),
   weib_shape0_5        = c(Hazard_Weib_C(shape.C.T.WEIB = 0.5,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Weib_T(shape.C.T.WEIB = 0.5,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   weib_shape1_5        = c(Hazard_Weib_C(shape.C.T.WEIB = 1.5,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Weib_T(shape.C.T.WEIB = 1.5,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   gomp_shape0_2        = c(Hazard_Gomp_C(shape.C.T.GOMP = 0.2,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Gomp_T(shape.C.T.GOMP = 0.2,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var)),
   gomp_shape_Minus_0_2 = c(Hazard_Gomp_C(shape.C.T.GOMP = -0.2,
                                          time_points = time_points,
                                          median.control = median.control),
                            Hazard_Gomp_T(shape.C.T.GOMP = -0.2,
                                          time_points = time_points,
                                          median.control = median.control,
                                          HR = HR, HR.var = HR.var))
) %>%
   pivot_longer(cols = c("exp", "weib_shape0_5", "weib_shape1_5",
                         "gomp_shape0_2", "gomp_shape_Minus_0_2"),
                names_to  = "dist",
                values_to = "hazard") %>%
   mutate(
      dist = factor(dist,
                    levels = c("exp", "weib_shape0_5", "weib_shape1_5",
                               "gomp_shape0_2", "gomp_shape_Minus_0_2"),
                    labels = c("Exponential",
                               "Weibull\n(shape = 0.5)", "Weibull\n(shape = 1.5)",
                               "Gompertz\n(shape = 0.2)", "Gompertz\n(shape = -0.2)"))
   )

# Plots
df_plot %>%
   filter(
      dist %in% c("Exponential", "Weibull\n(shape = 1.5)",
                  "Gompertz\n(shape = 0.2)"
      ) & grp == "Treatment"
   ) %>%
   ggplot(aes(x=t, y=hazard, col = dist)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Hazard function") +
   scale_color_manual("Failure time\ndistribution",
                      values = c("black", "red", "blue")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> hazard_increasing


df_plot %>%
   filter(
      dist %in% c("Exponential", "Weibull\n(shape = 0.5)",
                  "Gompertz\n(shape = -0.2)"
      ) & grp == "Treatment"
   ) %>%
   ggplot(aes(x=t, y=hazard, col = dist)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Hazard function") +
   scale_color_manual("", values = c("black", "orangered4", "royalblue4")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> hazard_decreasing

ggarrange(hazard_increasing, hazard_decreasing) -> AppendixFigure1_hazards
AppendixFigure1_hazards


#----------- APPENDIX: Survival function of late treatment effect --------------
# Functions:
Survival_Treatment <- function(x, l_c, l_t, c){
   # x:   time points
   # l_c: exponential parameter of control group (and treatment group until c)
   # l_t: exponential parameter of treatment group from c onwoards
   # c:   starting point of treatment effect in treatment group
   return(
      ifelse(x <= c,
             exp(-l_c*x),
             exp(-l_c*c) * exp(-l_t*(x-c)))
   )
}

# Assumed parameters
HR <- 0.7
HR.var <- 1
median.control <- 12
start_time <- 1/3*median.control
time_points <- seq(1,
                   50,
                   by = 0.2)
df_plot_Late_treatment <- data.frame(
   t            = rep(time_points, times=2),
   grp          = c(rep("Control: exponential failure times",
                        length(time_points)),
                    rep("Treatment: piece-wise exponential failure times",
                        length(time_points))),
   non.prop.exp = c(1-pexp(q=time_points, rate = log(2)/median.control),
                    Survival_Treatment(x=time_points,
                                       l_c = log(2)/median.control,
                                       l_t = (HR*HR.var)*log(2)/median.control,
                                       c=start_time))
) %>%
   mutate(non.prop.exp_100 = non.prop.exp*100)

# Data frame for plot
df_plot_Late_treatment %>%
   ggplot(aes(x=t, y=non.prop.exp_100, col = grp)) +
   geom_line(size=0.4) +
   labs(x = "Time (in months)", y = "Survival function") +
   scale_y_continuous(breaks = seq(0,100, by=20)) +
   scale_color_manual("Group", values = c("black", "red", "blue")) +
   theme_bw() +
   theme(
      text = element_text(size = 8,
                          family = "serif" # TT Times New Roman
                          #family = "sans" # TT Arial
      ),
      legend.position = "bottom",
      axis.ticks = element_line(size=0.4),
      axis.ticks.length = unit(0.05, "cm")
   ) -> AppendixFigure2_Survival_late_treatment
AppendixFigure2_Survival_late_treatment
