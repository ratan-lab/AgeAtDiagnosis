## Last updated: July 24, 2026 
## Writer: Jisu Shin


## This code is for Table 1. 

set.seed(123)
library(tidyverse)

## load data 

# cli = read_csv(file.path("//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/Age_heejin_Nate_v4_07102026.csv"), col_types = cols(.default = col_character()))
cli = read_csv(file.path("//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/clinical_final.csv"), col_types = cols(.default = col_character())) # updated version with Katie's cohort 

cli %>% dplyr::count(LGL_Type) %>% mutate(per = n/sum(n)*100)
# subset only T 
cli = cli %>% filter(LGL_Type == "T")

  # ------------ total number of patients  -------------------

  N = nrow(cli)
  N 

  # ------------ Age at diagnosis -------------------
  age_dx  <- cli$Age_dx %>% as.numeric()


  mean(age_dx) 
  sd(age_dx) 
  median(age_dx) 
  min(age_dx) 
  max(age_dx) 
  summary(age_dx)

# ------------ Sex, LGL_subtype,RA,STAT3_status ------------------

  cli %>% dplyr::count(Sex) %>% mutate(per = n/sum(n)*100)
  
  cli %>% filter(!is.na(RA_status)) %>% dplyr::count(RA_status) %>% mutate(per = n/sum(n)*100)
  cli %>% filter(!is.na(STAT3_status))  %>% dplyr::count(STAT3_status) %>% mutate(per = n/sum(n)*100, total_n = sum(n))
  cli %>% filter(!is.na(STAT3_status))  %>% dplyr::count(STAT3_mutations) %>% mutate(per = n/sum(n)*100)
  cli %>% filter(!is.na(STAT3_status) & STAT3_status == "MT")  %>% dplyr::count(STAT3_mutations) %>% mutate(per = n/sum(n)*100)


# ------------ CBC parameters  ------------------
    cbc_par = colnames(cli)[13:23]
    
    cli_filt = cli %>% filter(cbc_date_flag == "pass" &  tx_flag =="pass")
    df = NULL 
    for (c in cbc_par) {
      var = cli_filt %>% pull(c) %>% as.numeric() 
      n_avail = sum(!is.na(var))
      med = median(var, na.rm =T) 
      q1 = quantile(var, 0.25, na.rm = TRUE)
      q3 = quantile(var, 0.75, na.rm = TRUE)
      a = c(c, n_avail, med, q1, q3) 
      df = rbind(df, a)
    }

    colnames(df) = c("CBC", "N", "Median", "Q1", "Q3")

# ------------ Ever treatment status  ------------------

  cli %>% filter(!is.na(ever_tx)) %>%  dplyr::count(ever_tx) %>% mutate(per = n/sum(n)*100) 

    # which treatment ? 
    cli %>% filter(ever_tx == "Y" & !is.na(Tx1)) %>% dplyr::count(Tx1) %>% mutate(perc = n/sum(n)*100) %>% arrange(-perc) 

    # ------------ Tx1 response rate  ------------------

    # available Tx response 
    cli %>% filter(ever_tx == "Y") %>% dplyr::count(Tx1_response) 

    cli$tx_res = cli %>% mutate(tx_res = ifelse(Tx1_response %in% c("UNK", "unk", "Unclear", "Unk", 0, "NE"), NA, ifelse(Tx1_response == "cr", "CR", 
      ifelse(str_detect(Tx1_response,"intolerance"), "Intolerance", ifelse(Tx1_response == "CR/CMR", "CR", Tx1_response))))) %>% 
      pull(tx_res)

    cli %>% filter(ever_tx == "Y" & !is.na(tx_res)) %>% dplyr::count(tx_res) %>% mutate(perc = n/sum(n)*100, total_n = sum(n)) %>% arrange(-perc)
    
    # cli %>% filter(ever_tx == "Y" & !is.na(Tx1)) %>% dplyr::count(Tx1, tx_res) %>% 
    #   group_by(Tx1) %>% mutate(perc=n/sum(n)*100) %>% 
    #   arrange(Tx1, -perc)

# ------------ Deceased information (survival)  ------------------

    cli %>%dplyr::count(Deceased_indicator) %>% mutate(perc = n/sum(n)*100)


#
#  genomic subsets 
    n_wes = sub(!is.na(cli$patient_id) & cli$patient_id != "") 