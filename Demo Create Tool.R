demo_create <- function(rag){
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
  
  rag_alt <- rag
  rag_alt <- rag_alt[,-c(1,2)]
  rag_alt$'Patient.Characteristic...FGFR.Status.' <- 
    case_when(
      rag_alt$Patient.Characteristic...FGFR.Status. == 'target' ~ 'FGFR +',
      rag_alt$Patient.Characteristic...FGFR.Status. == 'wt' ~ 'FGFR -'
    )
  
  pan_tumor <- rag[which(rag$Patient.Characteristic..Pan.tumor.pt..),]
  specific_tumor <- rag[which(rag$Patient.Characteristic..Tumor.specific.pt.),]
  
  pan_tumor$'Patient.Characteristic...FGFR.Status.' <- 
    case_when(
      pan_tumor$Patient.Characteristic...FGFR.Status. == 'target' ~ 'FGFR +',
      pan_tumor$Patient.Characteristic...FGFR.Status. == 'wt' ~ 'FGFR -'
    )
  
  specific_tumor$'Patient.Characteristic...FGFR.Status.' <- 
    case_when(
      specific_tumor$Patient.Characteristic...FGFR.Status. == 'target' ~ 'FGFR +',
      specific_tumor$Patient.Characteristic...FGFR.Status. == 'wt' ~ 'FGFR -'
    )
  
  seperate_outcomes <- list(rag_alt, pan_tumor, specific_tumor)
  
  Typer <- function(df) {
    n <- ncol(df)
    
    type <- rep(NA, n)
    for(k in 1:n) {
      data_temp <- df[, k]
      if (is.numeric(data_temp[1])) {
        if (length(unique(data_temp)) < 10) {
          type[k] <- 'cat'
        }else {
          type[k] <- 'num'
        }
      }else {
        if (length(unique(data_temp)) < 10) {
          type[k] <- 'cat'
        }else {
          type[k] <- 'nr'
        }
      }
    }
    
    return(type)
  }
  
  types <- Typer(rag_alt)
  
  Tables_Create <- function(type, df, var, strat) {
    if (type == 'cat') {
      var <- df[, var]
      strat <- df[, strat]
      
      var[which(is.na(var))] <- 'NA'
      var[which(var == '')] <- 'NA'
      
      table_temp <- as.data.frame(table(var, strat))
      table_temp <- spread(table_temp, key = strat, value = Freq)
      
      n <- ncol(table_temp) - 1
      for(k in 2:(n+1)) {
        table_temp[,k] <-
          paste0(table_temp[,k],
                 ' (',
                 round(table_temp[,k]/sum(table_temp[,k]),3)*100,
                 '%',
                 ')')
      }
      if('NA' %in% table_temp$var){
        table_temp <- rbind(table_temp[which(table_temp$var != 'NA'),],
                            table_temp[which(table_temp$var == 'NA'),])
      }
      
      names(table_temp)[1] <- ''
      
    }else if(type == 'num'){
      table_temp <- as.data.frame(
        df %>% 
          group_by(.data[[strat]]) %>% 
          summarise(mean = mean(.data[[var]], na.rm = TRUE),
                    sd = sd(.data[[var]], na.rm = TRUE),
                    median = median(.data[[var]], na.rm = TRUE),
                    lower = quantile(.data[[var]], c(.25), na.rm = TRUE),
                    upper = quantile(.data[[var]], c(.75), na.rm = TRUE),
                    min = min(.data[[var]], na.rm = TRUE),
                    max = max(.data[[var]], na.rm = TRUE),
                    na = sum(is.na(.data[[var]])),
                    n = n()) %>% 
          mutate(`Mean (sd)` = paste0(round(mean, 2), ' (', round(sd, 2), ')')) %>% 
          mutate(`Median [IQR]` = paste0(round(median,2), ' [', round(lower, 2), ', ', round(upper, 2), ']')) %>% 
          mutate(`Range` = paste0(round(min, 2), ' - ', round(max, 2))) %>% 
          mutate(`Unknown; n (%)` = paste0(na,' (',
                                           round(na/n,3)*100,'%',')')) %>% 
          select( .data[[strat]], `Mean (sd)`, `Median [IQR]`, `Range`, `Unknown; n (%)`) %>% 
          t()
      ) 
      
      names(table_temp) <- table_temp[1,]
      table_temp <- table_temp[2:nrow(table_temp),]
    }
    
    return(table_temp)
  }
  
  count = 0
  for(k in c(1:5, 7:ncol(rag_alt))) {
    if (types[k] != 'nr') {
      count <- count + 1
      
      var <- names(rag_alt)[k]
      names_clean <- names(janitor::clean_names(rag_alt))[k]
      strat <- 'Patient.Characteristic...FGFR.Status.'
      
      for(g in 1:length(seperate_outcomes)){
        if (g == 1){
          demo_table <- Tables_Create(types[k], seperate_outcomes[[g]], var, strat)
        }else{
          if(types[k] == 'cat'){
            tab_temp <-  Tables_Create(types[k], seperate_outcomes[[g]], var, strat)
            demo_table <- merge(demo_table, tab_temp, by = c("", ""), all = TRUE)
            names(demo_table)[1] <- ''
          }else{
            tab_temp <-  Tables_Create(types[k], seperate_outcomes[[g]], var, strat)
            demo_table <- cbind(demo_table, tab_temp)
          }
        }
      }
      
      demo_table <-apply(demo_table, c(1,2), function(x) if(is.na(x)) {'0 (0%)'}else{x} ) 
      
      if(ncol(demo_table) == 7){
        kab_table <- knitr::kable(demo_table, caption = 
                                    paste0('<strong>', count, '. ', names_clean, '</strong>'),
                                  col.names = c('', rep(c('FGFR -', 'FGFR +'), 3))) %>%
          kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
          column_spec(1, bold = T, border_right = T, width = "20em") %>%
          column_spec(2:6, border_right = T, width = "20em") %>%
          column_spec(7, width = "20em") %>%
          add_header_above(c(" ", "Pan and Specific" = 2, "Pan" = 2, "Specific" = 2))
      }else{
        kab_table <- knitr::kable(demo_table, caption = 
                                    paste0('<strong>', count, '. ', names_clean, '</strong>'),
                                  col.names = c(rep(c('FGFR -', 'FGFR +'), 3))) %>%
          kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>% 
          column_spec(1, bold = T, border_right = T, width = "10em") %>%
          column_spec(2:6, border_right = T, width = "20em") %>%
          column_spec(7, width = "20em") %>%
          add_header_above(c(" ", "UDM Patients + CG Database" = 2, "UDM Patients" = 2, "CG Database" = 2))
      }
      
      print(kab_table)
      cat('\n')
    }
  }
}

setwd('C:\\Users\\matthew.baranoff\\Working\\RAGNAR\\Tools')
save(demo_create, file = 'demo_create.RData')