# 
# raganar Single Tumor Type Analysis
# Author: Matthew Baranoff
# Date: 3-14-22
# Access Date: 3-9-22
# 
# QC Author:
# QC Date:
# QC Comment:
# 
# 

## Load Needed Packages

## Analysis Function
RagnarOutcomes <- function(rag, name) {
  ## Load Needed Packages
  require(ggplot2)
  require(survminer)
  require(dplyr)
  require(tidyr)
  require(tidyverse)
  require(readxl)
  require(survival)
  require(ggpmisc)
  require(ggforce)
  require(broom)
  
  analysis <- list()
  
  #######################################################
  #                                                     #               
  #                 Define Functions                    #
  #                                                     #
  #######################################################
  
  ## Function used to set any value with 0.5 to 0
  SetZero <- function(df) {
    df[] <- apply(df, 
                  MARGIN = c(1, 2),
                  FUN = function(x) {if (x == 0.5) {0}else {x}})
    return(df)
  }
  
  ## Takes two vectors df_logical and df_numeric and return each value that has a true as the df_numeric otherwise return 9999
  ChooseTrue <- function(df_logical, df_numeric) {
    for(k in 1:ncol(df_logical)) {
      df_numeric[,k] <-mapply(function(x, y) {if (x) {y}else {9999}},  # defining 9999 as NA for logic later
                              df_logical[,k], 
                              df_numeric[,k])
    }
    return(df_numeric)
  }
  
  #######################################################
  #                                                     #               
  #            0 Loading and Cleaning                   #
  #                                                     #
  #######################################################
  ## Analysis Function
  
  ## Dataset to be worked on
  rag$'Patient.Characteristic...FGFR.Status.' <- 
    case_when(
      rag$Patient.Characteristic...FGFR.Status. == 'target' ~ 'FGFR +',
      rag$Patient.Characteristic...FGFR.Status. == 'wt' ~ 'FGFR -'
    )
  
  rag <- rag[which(rag$Patient.Characteristic...Disease. == 'breast'),]
  ## Define outcomes
  outcomes <- data.frame(last_day_clincal = rag$Patient.Outcome..CENSORING.VARIABLE..Last.day.of.clinical.activity..,
                         selective_FGFR = rag$Patient.Outcome..CENSORING.VARIABLE..Selective.FGFR.Inhibitor.use..,
                         study_drug = rag$Patient.Outcome..CENSORING.VARIABLE..Clinical.Study.Drug.Use..,
                         end_of_data = rag$Patient.Outcome..CENSORING.VARIABLE..End.of.Data..December.31..2020...,
                         death = rag$Patient.Outcome..OUTCOME..Death..)
  
  outcomes_days <- data.frame(last_day_clincal = rag$person.time.for.PO3,
                              selective_FGFR = rag$person.time.for.PO4,
                              study_drug = rag$person.time.for.PO5,
                              end_of_data = rag$person.time.for.PO6,
                              death = rag$person.time.for.PO7)
  
  #check if any censoring occurs after endofdata -- 0 people
  nstuff <- nrow(outcomes_days[which(outcomes_days$end_of_data < (outcomes_days$death | 
                                                                    outcomes_days$selective_FGFR | 
                                                                    outcomes_days$end_of_data | 
                                                                    outcomes_days$end_of_data)),])
  if(nstuff > 0){
    warning('Censoring occurs after end of data')
  }else{
    warning('Censoring does not occur after end of data')
  }
  
  ## Define delaying test
  delay <- data.frame(first_visit = rag$Patient.Outcome..IF.DELAYED.ENTRY..FGFR.Test.Result....or....,
                      second_visit = rag$Patient.Outcome..DELAYED.ENTRY.CRITERIA..Second.Visit..)
  
  delay_days <- data.frame(first_test = rag$person.time.for.PO1,
                           second_visit = rag$person.time.for.PO2)
  
  #check for visit after death -- 0 people
  nstuff <- nrow(outcomes_days[which(outcomes_days$death < delay_days$second_visit),])
  if(nstuff > 0){
    warning('Second visit occurs after death')
  } else {
    warning('Second visit does not occurs after death')
  }
  
  ## define var type for demographics
  
  ## Apply logic to outcomes and delays 0.5 values <- 0 + choosing true values
  outcomes_days <- SetZero(outcomes_days)
  delay_days <- SetZero(delay_days)
  
  outcomes_days <- ChooseTrue(outcomes, outcomes_days)
  delay_days <- ChooseTrue(delay, delay_days)
  
  ## Death at 0 <- death + 16
  outcomes_days[which(outcomes_days$death == 0), 'death'] <- 16
  
  #check where death occurs before last clinical activity - 545 people
  nstuff <- nrow(outcomes_days[which(outcomes_days$last_day_clincal > outcomes_days$death),])
  warning(paste('Death occurs before the last day of clincal activity in', nstuff, 'cases', sep = ' '))
  
  ## Logic of last visit
  outcomes_days <- outcomes_days %>%
    mutate(last_day_clincal = 
             ifelse(last_day_clincal < death & death < end_of_data, death,
                    ifelse(last_day_clincal < end_of_data & end_of_data < death, end_of_data,
                           9999))
    )
  ## Apply censoring logic
  outcomes_days$censored_date <- 
    apply(outcomes_days[,1:(ncol(outcomes_days)-1)], 1, min)
  
  ## If death after censord value choose censored date, 1 is defined as death, 0 for censored, no survivors
  outcomes_days$outcome_defined <-
    ifelse(outcomes_days$death > outcomes_days$censored_date | outcomes_days$death == 9999, 0, 1)
  
  ## Define number of days to censor
  outcomes_days$days <-
    ifelse(outcomes_days$outcome_defined == 0, outcomes_days$censored_date, outcomes_days$death)
  
  ## Applying delayed entry logic to those with delayed entry
  delay_days$entry <- ifelse((delay$first_visit  & !rag$Patient.Characteristic..BEFORE.OR.ON.CED..FGFR.Test.Result....or....) |
                               (delay$second_visit & !rag$Patient.Characteristic..BEFORE.CED..Second.Visit.),
                             TRUE, FALSE)
  
  delay_days$first_test <- ifelse((delay$first_visit  & !rag$Patient.Characteristic..BEFORE.OR.ON.CED..FGFR.Test.Result....or....), 
                                  delay_days$first_test, 9999)
  
  delay_days$second_visit <- ifelse((delay$second_visit  & !rag$Patient.Characteristic..BEFORE.CED..Second.Visit.), 
                                    delay_days$second_visit, 9999)
  
  ## Make NA's negative
  delay_days[] <- apply(delay_days, 
                        MARGIN = c(1, 2),
                        FUN = function(x) {if (x == 9999) {-9999}else {x}})
  
  ## Defined dayes to delay entry
  delay_days$day <- apply(delay_days[,1:2], 1, FUN = max)
  delay_days[which(!delay_days$entry), 'day'] <- 0
  
  ## Removed delayed entry from days to event
  outcomes_days$days_delay <- delay_days$day
  outcomes_days$is_delayed <- ifelse(delay_days$day == 0, 0, 1)
  
  nstuff <- nrow(outcomes_days[which(outcomes_days$days_delay > rag$person.time.for.PO7),])
  warning(paste(nstuff, 'observation(s) have a test after death', sep = ' '))
  
  #######################################################
  #                                                     #               
  #        Cox PH models Without Delayed Entry          #
  #                                                     #
  #######################################################
  
  ## First model without any adjustments
  cox_unadjusted_data <- cbind(FGFR = rag$Patient.Characteristic...FGFR.Status., outcomes_days)
  
  ## Define the referant group
  cox_unadjusted_data$FGFR = factor(cox_unadjusted_data$FGFR, 
                                    levels = c("FGFR -","FGFR +"))
  
  cox_unadjusted <- coxph(Surv(days, outcome_defined) ~ FGFR + strata(rag$Patient.Characteristic...Disease.), data = cox_unadjusted_data)
  
  analysis[['cox_unadjusted']] <- cox_unadjusted 
  
  ## Minimally adjusted model
  cox_min_adjusted_data <- cbind(FGFR = rag$Patient.Characteristic...FGFR.Status., 
                                 outcomes_days,
                                 age_meta = rag$Patient.Characteristic..DEMO..2c..Age.at.metastatic.diagnosis..tumor.specific...pan.tumor...,
                                 sex = rag$Patient.Characteristic..DEMO..4..S.,
                                 baitset = rag$Patient.Characteristic..Baitset...Recat..,
                                 date_of_meta = rag$Patient.Characteristic..FINAL..Year.of.Adv.Met.Diagnosis.,
                                 cancer_type = rag$Patient.Characteristic...Disease.)
  
  ## Replace cancers with < 10 obs as other
  cancer_num <- as.data.frame(table(cox_min_adjusted_data$cancer_type))
  for(k in 1:nrow(cox_min_adjusted_data)) {
    temp_cancer <- cox_min_adjusted_data$cancer_type[k]
    temp_cancer_num <- cancer_num[which(cancer_num$Var1 == temp_cancer), 'Freq']
    if (temp_cancer_num < 10) {
      cox_min_adjusted_data$cancer_type[k] <- 'Other'
    }
  }
  
  ## Choose referent type as wt
  cox_min_adjusted_data$FGFR = factor(cox_min_adjusted_data$FGFR, 
                                      levels = c("FGFR -","FGFR +"))
  
  ## Remove sex from cancers that inlcude only one sex, same for cancers
  if (length(unique(cox_min_adjusted_data$sex)) == 1) {
    cox_min_adjusted_data <- cox_min_adjusted_data %>% 
      select(-'sex')
  }
  if (length(unique(cox_min_adjusted_data$cancer_type)) == 1) {
    cox_min_adjusted_data <- cox_min_adjusted_data %>% 
      select(-'cancer_type')
  }
  
  ## Set date_of_meta to factor level variable
  cox_min_adjusted_data <- cox_min_adjusted_data %>% 
    mutate(date_of_meta = factor(date_of_meta))
  
  ## Run Min Adjusted model
  cox_min_adjusted <- coxph(Surv(cox_min_adjusted_data$days, cox_min_adjusted_data$outcome_defined) ~ . + strata(cox_min_adjusted_data$cancer_type), data =  
                              cox_min_adjusted_data[,c(1,11:ncol(cox_min_adjusted_data))])
  analysis[['cox_min_adjusted']] <- cox_min_adjusted
  
  ## fully adjusted model
  cox_fully_adjusted_data <- cbind(FGFR = rag$Patient.Characteristic...FGFR.Status., 
                                   outcomes_days,
                                   age_meta = rag$Patient.Characteristic..DEMO..2c..Age.at.metastatic.diagnosis..tumor.specific...pan.tumor...,
                                   sex = rag$Patient.Characteristic..DEMO..4..S.,
                                   baitset = rag$Patient.Characteristic..Baitset...Recat..,
                                   date_of_meta = rag$Patient.Characteristic..FINAL..Year.of.Adv.Met.Diagnosis.,
                                   immu_flag = rag$Patient.Characteristic..For.DEMO..15..IT..,
                                   chemo_flag = rag$Patient.Characteristic..For.DEMO..15..CT..,
                                   horm_flag = rag$Patient.Characteristic..For.DEMO..15..HT..,
                                   cci = rag$Patient.Characteristic..DEMO..23..CCI.Score.Measure..,
                                   group_stage = rag$Patient.Characteristic..DEMO..10..Group.Stage..,
                                   smoking = rag$Patient.Characteristic..DEMO..6..Smoking.History..,
                                   cancer_type = rag$Patient.Characteristic...Disease.)
  
  ## Replce cancers with < 10 obs as other
  cancer_num <- as.data.frame(table(cox_fully_adjusted_data$cancer_type))
  
  for(k in 1:nrow(cox_fully_adjusted_data)) {
    temp_cancer <- cox_fully_adjusted_data$cancer_type[k]
    temp_cancer_num <- cancer_num[which(cancer_num$Var1 == temp_cancer), 'Freq']
    
    if (temp_cancer_num < 10) {
      cox_fully_adjusted_data$cancer_type[k] <- 'Other'
    }
  }
  
  ## Choose referent type as wt
  cox_fully_adjusted_data$FGFR = factor(cox_fully_adjusted_data$FGFR, 
                                        levels = c("FGFR -","FGFR +"))
  
  ## Remove sex from cancers that inlcude only one sex, same for cancers
  if (length(unique(cox_fully_adjusted_data$sex)) == 1) {
    cox_fully_adjusted_data <- cox_fully_adjusted_data %>% 
      select(-'sex')
  }
  if (length(unique(cox_fully_adjusted_data$cancer_type)) == 1) {
    cox_fully_adjusted_data <- cox_fully_adjusted_data %>% 
      select(-'cancer_type')
  }
  
  ## Set date_of_meta to factor level variable
  cox_fully_adjusted_data <- cox_fully_adjusted_data %>% 
    mutate(date_of_meta = factor(date_of_meta)) %>% 
    mutate(group_stage = factor(group_stage))
  
  ## Run Min Adjusted model
  cox_fully_adjusted <- coxph(Surv(cox_fully_adjusted_data$days, cox_fully_adjusted_data$outcome_defined) ~ . + strata(cox_fully_adjusted_data$cancer_type), data =  
                                cox_fully_adjusted_data[,c(1,11:ncol(cox_fully_adjusted_data))])
  analysis[['cox_fully_adjusted']] <- cox_fully_adjusted 
  
  #######################################################
  #                                                     #               
  #        Cox PH models With Delayed Entry             #
  #                                                     #
  #######################################################
  
  ## First model without any adjustments
  cox_unadjusted_data_delay <- cox_unadjusted_data[which(cox_unadjusted_data$days -  cox_unadjusted_data$days_delay > 0),]
  
  cox_unadjusted_delayed <- coxph(Surv(time = days_delay, 
                                       time2 = days, 
                                       event = outcome_defined) ~ FGFR + strata(rag$Patient.Characteristic...Disease.[which(cox_unadjusted_data$days -  cox_unadjusted_data$days_delay > 0)]), 
                                  data = cox_unadjusted_data_delay)
  analysis[['cox_unadjusted_delayed']] <- cox_unadjusted_delayed 
  
  ## Minimally adjusted model
  ## Remove sex from cancers that inlcude only one sex
  cox_min_adjusted_data_delay <- cox_min_adjusted_data[which(cox_min_adjusted_data$days -  cox_min_adjusted_data$days_delay > 0),]
  cox_min_adjusted_data_delay <- na.omit(cox_min_adjusted_data_delay)
  
  cox_min_adjusted_delayed <-  coxph(Surv(time = cox_min_adjusted_data_delay$days_delay, 
                                          time2 = cox_min_adjusted_data_delay$days, 
                                          event = cox_min_adjusted_data_delay$outcome_defined) ~ . + strata(cox_min_adjusted_data_delay$cancer_type), 
                                     data = cox_min_adjusted_data_delay[,c(1,11:(ncol(cox_min_adjusted_data_delay)-1))])
  
  analysis[['cox_min_adjusted_delayed']] <- cox_min_adjusted_delayed 
  
  ## Fully adjusted model
  ## Remove sex from cancers that inlcude only one sex
  cox_fully_adjusted_data_delay <- cox_fully_adjusted_data[which(cox_fully_adjusted_data$days -  cox_fully_adjusted_data$days_delay > 0),]
  cox_fully_adjusted_data_delay <- na.omit(cox_fully_adjusted_data_delay)
  
  cox_fully_adjusted_delayed <-  coxph(Surv(time = cox_fully_adjusted_data_delay$days_delay, 
                                            time2 = cox_fully_adjusted_data_delay$days, 
                                            event = cox_fully_adjusted_data_delay$outcome_defined) ~ .  + strata(cox_fully_adjusted_data_delay$cancer_type), 
                                       data = cox_fully_adjusted_data_delay[,c(1,12:(ncol(cox_fully_adjusted_data_delay)-1))])
  analysis[['cox_fully_adjusted_delayed']] <- cox_fully_adjusted_delayed
  
  #######################################################
  #                                                     #               
  #                     KM plots                        #
  #                                                     #
  #######################################################
  
  # KMcreate <- function(rag, time_var, name_var, plot_name, delay = FALSE) {
  #   
  #   if(!delay){
  #     km_data  <- summary(survfit(Surv(outcomes_days$days, outcomes_days$outcome_defined) ~ 
  #                                   rag$Patient.Characteristic...FGFR.Status.))
  #   }
  #   
  #   km <- data.frame(time = km_data$time, survival = km_data$surv)
  #   km$type <- rep(NA, nrow(km))
  #   
  #   start <- 'FGFR-' 
  #   km$type[1] <- start
  #   for(k in 2:nrow(km)){
  #     if(km$time[k] < km$time[k-1]){
  #       start <- 'FGFR+'
  #     }
  #     km$type[k] <- start
  #   }
  #   
  #   cen <- outcomes_days[which(outcomes_days$outcome_defined == 0),]
  #   cen <- cbind(cen, FGFR = rag$Patient.Characteristic...FGFR.Status.[which(outcomes_days$outcome_defined == 0)])
  #   cen$risk <- rep(NA, nrow(cen))
  #   
  #   for(k in 1:nrow(cen)){
  #     max_t <- max(km$time[which(km$time < cen$days[k] & km$type)])
  #     cen$risk[k] <- km$survival[km$time == max_t]
  #   }
  #   
  #   ## Define the limits of the plots
  #   limit = max(ceiling(km_data$years))
  #   
  #   ## Create KM plot
  #   KMmain <- ggplot(km_data, aes(x = years, y = surv_prop, group = group)) +
  #     geom_step(aes(color = group), alpha = 0.25, size = 1.5)  +
  #     geom_point(data = censored_km_data, aes(x = years, y = surv_prop, color = group),
  #                shape=3, size = 2) +
  #     theme_minimal() +
  #     xlab('Years') +
  #     scale_y_continuous(labels = scales::percent, limits=c(0,1)) +
  #     ylab("Survival Probability") +
  #     scale_color_manual(values=c('red','blue'), name = "FGFR Type") +
  #     theme(legend.position="top", text = element_text(size = 17)) +
  #     scale_x_continuous(n.breaks = 10, limits = c(0,limit)) +
  #     ggtitle(plot_name) +
  #     theme(text=element_text(face="bold", size=12))
  #   
  #   ## Create data table
  #   posCenval <- seq(0, limit - 1)
  #   negCenval <- seq(0, limit - 1)
  #   
  #   ## Find the number censored for each year cutoff by FGFR type
  #   for(k in 1:limit) {
  #     posCenval[k] <- length(which(km_data$group == 'FGFR +' & km_data$years < k & 
  #                                    km_data$outcome == 0 &  km_data$years > 0))
  #     negCenval[k] <- length(which(km_data$group == 'FGFR -' & km_data$years < k & 
  #                                    km_data$outcome == 0 & km_data$years > 0))
  #   }
  #   
  #   posval <- seq(0, limit - 1)
  #   negval <- seq(0, limit - 1)
  #   
  #   ## Find the number at risk for each year by FGFR type 
  #   for(k in 1:limit) {
  #     posval[k] <- length(which(km_data$group == 'FGFR +' & km_data$years > k & km_data$years > 0))
  #     negval[k] <- length(which(km_data$group == 'FGFR -' & km_data$years > k & km_data$years > 0))
  #   }
  #   
  #   ## Create a table of the at risk and by FGFR type
  #   km_dataTable <- data.frame(values = c(groups$Freq[2], 
  #                                         posval,
  #                                         groups$Freq[1], 
  #                                         negval),
  #                              x = rep(seq(0, limit),2), 
  #                              y = c(rep('FGFR+ at Risk \n (Censored)', limit + 1), rep('FGFR- at Risk \n (Censored)', limit + 1)))
  #   
  #   ## Create 'death (censored)' string
  #   km_dataTable$values[which(km_dataTable$y == 'FGFR+ at Risk \n (Censored)')][2:(limit + 1)] <- paste0(
  #     km_dataTable$values[which(km_dataTable$y == 'FGFR+ at Risk \n (Censored)')][2:(limit + 1)], '\n',
  #     " (", posCenval[1:limit], ")"
  #   )
  #   
  #   km_dataTable$values[which(km_dataTable$y == 'FGFR- at Risk \n (Censored)')][2:(limit + 1)] <- paste0(
  #     km_dataTable$values[which(km_dataTable$y == 'FGFR- at Risk \n (Censored)')][2:(limit + 1)], '\n',
  #     " (", negCenval[1:limit], ")"
  #   )
  #   
  #   ## Draw KM table
  #   KMtable <- ggplot(km_dataTable, aes(x = y, y = x)) + theme_minimal() + 
  #     annotate("text", y = km_dataTable$x[1:(limit + 1)], x = km_dataTable$y[1:(limit + 1)], 
  #              label = km_dataTable$values[1:(limit + 1)], color = '#0000FF',
  #              size = 6) + 
  #     coord_flip() +
  #     annotate("text", y = km_dataTable$x[(limit + 2):(2*(limit + 1))], x = km_dataTable$y[(limit + 2):(2*(limit + 1))], 
  #              label = km_dataTable$values[(limit + 2):(2*(limit + 1))], color = '#FF0000',
  #              size = 6) + 
  #     coord_flip() + geom_hline(yintercept = seq(0.5, 8.5, by = 1)) +
  #     geom_vline(xintercept = 1.5) +
  #     theme(panel.grid  = element_blank(), text = element_text(size = 17),
  #           axis.text.x=element_blank()) +
  #     ylab('') + xlab("") + ylim(0,limit)
  #   
  #   ## Create and save KM table and KM plot
  #   return(ggarrange(KMmain, KMtable, heights = c(2, 1), 
  #                    ncol = 1, nrow = 2, align = 'v'))
  # }
  # 
  # analysis[[paste0("No_Delay", name)]] <- 
  #   KMcreate(rag, outcomes_days$days, none, 
  #            str_replace_all(paste0("No Delay ", name, Sys.Date()), "_", " "))
  # 
  # outcomes_days <- outcomes_days[which(outcomes_days$days > outcomes_days$days_delay),]
  # analysis[[paste0("With_Delay", name)]] <- 
  #   KMcreate(rag, outcomes_days$days, none, 
  #            str_replace_all(paste0("With  Delay ", name, Sys.Date()), "_", " "), delay = TRUE)
  
  #check where death occurs after endofdata - 0 people
  if(nrow(outcomes_days[which(outcomes_days$end_of_data < outcomes_days$death),]) > 0){
    warning('End of Data Occurs Before Death')
  }else{
    warning('End of Data Does not Occurs Before Death')
  }
  
  #histogram of delay days
  print(outcomes_days %>%
          ggplot(aes(x = days_delay, fill=rag$Patient.Characteristic...FGFR.Status.)) +
          geom_histogram( color="#E9ECEF", alpha=0.6, position = 'identity', bins=100) +
          scale_fill_manual(values=c("#69B3A2", "#404080")) +
          labs(fill=""))
  
  return(analysis)
}

setwd('C:\\Users\\matthew.baranoff\\Working\\RAGNAR\\Tools')
save(RagnarOutcomes, file = 'RagnarOutcome.RData')

