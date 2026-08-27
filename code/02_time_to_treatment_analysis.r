

## Date: July 29, 2026 
## writer: Jisu Shin 

## Time to treatment and survival analysis for T-LGLL. 


# ----------------- Load packages and data --------------------
library(data.table)
library(tidyverse)
library(lubridate)
library(ggpubr) 
library(haven)
library(cowplot) 
library(cocor) 
library(readxl)
library(survival)
library(survminer)
library(patchwork)

my_theme <- theme_minimal(base_size = 7) +
  theme(
    plot.title      = element_text(size = 8),
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 7),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    strip.text      = element_text(size = 8),
    text            = element_text(family = "Arial")
  )



fig.path ="//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/figures"


# cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/Age_heejin_Nate_v4_07102026.csv") 
cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/clinical_final.csv") 
cli = cli %>% filter(LGL_Type == "T")

# ----------------- 0. prep data   --------------------

cli_qced <- cli %>% 
    mutate(Age_dx = as.numeric(Age_dx),
    Sex = factor(Sex, levels = c("F", "M")),
    STAT3_status = factor(STAT3_status, levels = c("WT", "MT")),
    time_to_treatment_qc_pass = as.logical(time_to_treatment_qc_pass),
    death_used_for_censor     = as.logical(death_used_for_censor)
  )

nrow(cli_qced)
table(cli_qced$time_to_treatment_qc_pass)
sum(cli_qced$time_to_treatment_qc_pass & cli_qced$ever_tx == "Y",na.rm = TRUE)
# median time to tx1 initiation date 
cli_qced %>% filter(time_to_treatment_qc_pass == TRUE & ever_tx == "Y") %>% pull(time_to_event) %>% median()
# median time to event date! 
median(cli_qced$time_to_event[cli_qced$time_to_treatment_qc_pass], na.rm = TRUE)

# fix 

cli_qced %>% filter(time_to_treatment_qc_pass == TRUE) %>% pull(time_to_event) %>% range()

# cli_qced %>% filter(Data_source == "Katie") %>% filter(time_to_treatment_qc_pass == TRUE) %>% 
#     dplyr::select(RegID, STAT3_status, ever_tx, dx_date, last_fu_raw, Tx1_date_raw,time_to_event) %>% 
#     mutate(event_date = ifelse(ever_tx == "Y", Tx1_date_raw, last_fu_raw)) %>%
#     mutate(event_date = lubridate::mdy(event_date),dx_date = lubridate::mdy(dx_date))  %>% 
#     mutate(time_to_event_check = as.numeric(event_date - dx_date)/30.44) %>%
#     dplyr::select(-event_date) %>% arrange(STAT3_status, ever_tx) %>% as.data.frame

## ── Master analytic dataset ───────────────────────────────────────────────────
cox.master <- cli_qced |>
  filter(time_to_treatment_qc_pass == TRUE) |>
  mutate(
    ever_tx_num = case_when(
      ever_tx == "Y" ~ 1L,
      ever_tx == "N" ~ 0L,
      TRUE           ~ NA_integer_
    ),
    RA_status = factor(
      case_when(RA == "Y" ~ "Y", RA == "N" ~ "N", TRUE ~ NA_character_),
      levels = c("N", "Y")
    ),
    cr_status = case_when(
      ever_tx == "Y"                                                          ~ 1L,
      ever_tx == "N" & !is.na(death_used_for_censor) & death_used_for_censor ~ 2L,
      TRUE                                                                    ~ 0L
    )
  ) |>
  dplyr::select(time_to_event, ever_tx_num, cr_status,
                Age_dx, Sex, STAT3_status, RA_status) |>  # 2 RA_status = NA, 1 patient with cr_staus == 2 (death before tx) --> total three patients removed for the analysis 
  filter(!is.na(ever_tx_num), !is.na(Age_dx), !is.na(Sex),
         !is.na(STAT3_status), !is.na(RA_status))

nrow(cox.master) 
sum(cox.master$ever_tx_num == 1)
median(cox.master$time_to_event, na.rm = TRUE)

# ----------------- 1. cox models   --------------------

    # basic model ~ age + sex
        main.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex, data = cox.master)
        print(summary(main.cox))
        print(cox.zph(main.cox))

    # STAT3 model; ~ age + sex + STAT3 
        stat3.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + STAT3_status,
                        data = cox.master)
        print(summary(stat3.cox))
        print(cox.zph(stat3.cox))

        stat3.km <- survfit(
            Surv(time_to_event, ever_tx_num) ~ STAT3_status,
            data = cox.master
        )

        summary(stat3.km)$table
         
         # log rank pvalue 
            stat3.logrank <- survdiff(
                Surv(time_to_event, ever_tx_num) ~ STAT3_status,
                data = cox.master
            )


    # RA model 
        ra.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + RA_status,
                    data = cox.master)
        print(summary(ra.cox))

        ra.km <- survfit(
            Surv(time_to_event, ever_tx_num) ~ RA_status,
            data = cox.master
        )

        summary(ra.km)$table
        ra.logrank <- survdiff(Surv(time_to_event, ever_tx_num) ~ RA_status, data = cox.master)
        # cox.master %>% group_by(RA_status) %>% summarise(med = median(time_to_event))
        
    # full model 
        full.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + STAT3_status + RA_status,
                        data = cox.master)
        print(summary(full.cox))

    # age*STAT3 interacton test 
    stat3.int <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx * STAT3_status + Sex,
                    data = cox.master)
    print(anova(stat3.cox, stat3.int, test = "LRT"))

    ra.int <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx * RA_status + Sex,
                    data = cox.master)
    cat("\nLRT: Age × RA interaction\n")
    print(anova(ra.cox, ra.int, test = "LRT"))

    cat("\nCompeting-risk event counts:\n")
        print(table(cox.master$cr_status, dnn = "0=censored/1=treatment/2=death"))


# ----------------- section 2. time to treatment related. but see whether the CBC parameters can predict the treatment need  --------------------

#  regression: CBC parameters predict treatment need
    
    cli_filt = cli %>% rename(Age = Age_Dx, STAT3 = STAT3_status)  %>% 
                filter(cbc_date_flag == "pass" & tx_flag == "pass") # 239 patients left 

    cd_cols = grep("^CD", colnames(cli_filt))
    for (g in cd_cols) cli_filt[[g]] <- suppressWarnings(as.numeric(cli_filt[[g]]))

    cli.new <- cli_filt %>% 
        mutate(
            STAT3    = factor(STAT3, levels = c("WT", "MT")),
            new.kit  = Clinic_Kit == "Kit",
            new.fem  = Sex == "F",
            new.ra   = RA == "Y",
            new.ever = Ever_Tx == "Y"
        )

    tjur_r2 <- function(model) {
        f <- fitted(model)
        y <- model$y
        mean(f[y == 1]) - mean(f[y == 0])
    }

    cbc_candidate_vars <- c("HGB", "HCT", "MCV", "RBC", "WBC",
                            "RDW", "Platelets", "ANC", "ALC")

    dat_full <- cli.new |>
    dplyr::select(new.ever, Age, all_of(cbc_candidate_vars)) |>
    na.omit()

    cat(sprintf("\nLogistic N = %d (%d treated, %d untreated)\n",
                nrow(dat_full), sum(dat_full$new.ever), sum(!dat_full$new.ever)))

    m_sat  <- glm(new.ever ~ ., data = dat_full, family = binomial)
    m_base <- step(m_sat, trace = 0)
    # m_base <- step(m_sat, trace = 1)

    step_vars <- names(coef(m_base))[-1]
    cat("Stepwise-retained predictors:", paste(step_vars, collapse = ", "), "\n")
    cat(sprintf("Tjur R² (CBC model) = %.3f\n", tjur_r2(m_base)))

    ci_base <- suppressMessages(confint(m_base))
    cat("\nCoefficients (CBC model):\n")
    print(round(cbind(beta = coef(m_base), ci_base, p = coef(summary(m_base))[, 4]), 4))

    ## add STAT3 to stepwise model
    dat_stat3 <- cli.new |>
    dplyr::select(new.ever, STAT3, all_of(step_vars)) |>
    filter(!is.na(STAT3)) |>
    na.omit()

    cat(sprintf("With STAT3 restriction: N = %d\n", nrow(dat_stat3)))

    formula_base  <- as.formula(paste("new.ever ~", paste(step_vars, collapse = " + ")))
    formula_stat3 <- update(formula_base, . ~ . + STAT3)
    m_base2       <- glm(formula_base,  data = dat_stat3, family = binomial)
    m_stat3       <- glm(formula_stat3, data = dat_stat3, family = binomial)

    lrt     <- anova(m_base2, m_stat3, test = "LRT")
    lrt_p   <- lrt[2, "Pr(>Chi)"]
    lrt_chi2 <- lrt[2, "Deviance"]
    cat(sprintf("LRT χ²(1) = %.3f, p = %.4f\n", lrt_chi2, lrt_p))
    cat(sprintf("Tjur R² base = %.3f | +STAT3 = %.3f\n",
                tjur_r2(m_base2), tjur_r2(m_stat3)))

ci_stat3 <- suppressMessages(confint(m_stat3))
s        <- coef(summary(m_stat3))
cat("\nSTAT3 coefficient:\n")
print(round(cbind(beta = s["STAT3MT", 1], ci_stat3["STAT3MT", , drop = FALSE],
                  p    = s["STAT3MT", 4]), 4))

# Figure 2A. STAT3 MT vs WT survival curve 
    km.data <- cox.master |>
            mutate(age60_group = factor(ifelse(Age_dx > 60, "Age > 60", "Age <= 60"),levels = c("Age <= 60", "Age > 60")))
    
    km.stat3 = km.data %>% dplyr::select(time_to_event, ever_tx_num, STAT3_status) %>% na.omit()
    fit.stat3 = survfit(Surv(time_to_event, ever_tx_num) ~ STAT3_status, data = km.stat3)
    print(survdiff(Surv(time_to_event, ever_tx_num) ~ STAT3_status, data = km.stat3))
    
    # manual p 
        stat3.logrank <- survdiff(
        Surv(time_to_event, ever_tx_num) ~ STAT3_status,
        data = km.stat3
        )

        stat3.p <- pchisq(
        stat3.logrank$chisq,
        df = length(stat3.logrank$n) - 1,
        lower.tail = FALSE
        )

        stat3.p

    p.km.stat3 <- ggsurvplot(
        fit.stat3, data = km.stat3,
        palette     = c("black", "grey60"),
        legend.labs = levels(km.stat3$STAT3_status),
        xlab        = "Time from diagnosis (months)",
        ylab        = "Treatment-free probability",
        # pval        = TRUE, 
        pval = sprintf("p = %.4f", stat3.p),
        pval.size = 3,
        conf.int = TRUE, conf.int.alpha = 0.15,
        risk.table  = FALSE,
        ggtheme     = my_theme
    )$plot + theme_classic() + 
        theme(legend.position = c(0.8,0.7), 
            legend.title = element_blank(), 
            legend.text = element_text(size = 7),
            legend.key.height = unit(5, "mm"),
            legend.key.width  = unit(5, "mm"),

            legend.spacing.y = unit(0.5, "mm"),
            legend.box.spacing = unit(1, "mm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 1))

    ggsave(file.path(fig.path, "Figure2A.png"), p.km.stat3, width = 8, height = 8, units = "cm", dpi = 300)

    # Figure S5. Age x RA_status
            km.ra  <- km.data |> dplyr::select(time_to_event, ever_tx_num, RA_status) |> na.omit()
            fit.ra <- survfit(Surv(time_to_event, ever_tx_num) ~ RA_status, data = km.ra)
                
            p.km.ra <- ggsurvplot(
            fit.ra, data = km.ra,
            palette     = c("black", "grey60"),
            legend.title = "RA status",
            legend.labs = levels(km.ra$RA_status),
            xlab        = "Time from diagnosis (months)",
            ylab        = "Treatment-free probability",
            pval        = TRUE, 
            pval.size = 3,
            conf.int = TRUE, conf.int.alpha = 0.15,
            risk.table  = FALSE,
            ggtheme     = my_theme
            )$plot + 
            theme_classic() + 
            theme(legend.position = c(0.8,0.7),
            legend.title = element_text(size =7), 
                legend.text = element_text(size = 7),
                legend.key.height = unit(5, "mm"),
                legend.key.width  = unit(5, "mm"),

                legend.spacing.y = unit(0.5, "mm"),
                legend.box.spacing = unit(1, "mm"),
                legend.margin = margin(t = 0, r = 0, b = 0, l = 1))

        ggsave(file.path(fig.path, "FigureS5.png"), p.km.ra, width = 8, height = 8, units = "cm", dpi = 300)

# Figure S6. Age x STAT3 four group comparisons 

    # when age 60; 
        km.age_stat3 <- km.data |>
            mutate(age60_stat3 = factor(paste(age60_group, STAT3_status, sep = " / "), levels = c("Age <= 60 / WT", "Age <= 60 / MT", "Age > 60 / WT", "Age > 60 / MT"))) |>
            dplyr::select(time_to_event, ever_tx_num, age60_stat3) |> na.omit()

        fit.age_stat3 <- survfit(Surv(time_to_event, ever_tx_num) ~ age60_stat3, data = km.age_stat3)
        print(fit.age_stat3)

        p.km.age_stat3 <- ggsurvplot(
            fit.age_stat3, data = km.age_stat3,
            palette     = c("black", "#D55E00", "grey60", "#E69F00"),
            legend.labs = levels(km.age_stat3$age60_stat3),
            xlab        = "Time from diagnosis (months)",
            ylab        = "Treatment-free probability",
            pval        = TRUE, pval.size = 3,
            conf.int = TRUE, conf.int.alpha = 0.15,
            risk.table  = FALSE,
            ggtheme     = my_theme
            )$plot & theme_classic() + 
            theme(legend.position = c(0.8,0.7), 
                legend.title = element_blank(), 
                legend.text = element_text(size = 7),
                legend.key.height = unit(5, "mm"),
                legend.key.width  = unit(5, "mm"),

                legend.spacing.y = unit(0.5, "mm"),
                legend.box.spacing = unit(1, "mm"),
                legend.margin = margin(t = 0, r = 0, b = 0, l = 1))
        
        ggsave(file.path(fig.path, "FigureS6.png"), p.km.age_stat3, width = 8, height = 8, units = "cm", dpi = 300)

    # finding the cutpoint for the four groups 
            # age_range <- km.data |> pull(Age_dx) |> quantile(c(0.15, 0.85), na.rm = TRUE) # remove outliers 
            # candidate_cutoffs <- seq(ceiling(age_range[1]), floor(age_range[2]), by = 1)

            # scan_result <- map_dfr(candidate_cutoffs, function(cut) {
            # dat <- km.data |>
            #     mutate(
            #     age_group    = ifelse(Age_dx <= cut, paste0("Age<=", cut), paste0("Age>", cut)),
            #     age_stat3_4g = factor(paste(age_group, STAT3_status, sep = " / "))
            #     ) |>
            #     dplyr::select(time_to_event, ever_tx_num, age_stat3_4g) |>
            #     na.omit()

            # grp_n <- table(dat$age_stat3_4g)
            # # test at least 10 patients per group  
            
            # if (length(grp_n) < 4 || min(grp_n) < 10) {
            #     return(tibble(cutoff = cut, p = NA_real_, min_group_n = min(grp_n)))
            # }

            # sd <- survdiff(Surv(time_to_event, ever_tx_num) ~ age_stat3_4g, data = dat)
            # p  <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
            # tibble(cutoff = cut, p = p, min_group_n = min(grp_n))
            # })

    #         p = ggplot(scan_result, aes(x = cutoff, y = p)) +
    #             geom_line() + geom_point() +
    #             geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    #             labs(x = "Age cutoff", y = "4-group log-rank p-value") + 
    #             my_theme + 
    #             theme_classic()
            
    #         ggsave(file.path(fig.path, "test.png"), p, width = 8, height = 8, units = "cm", dpi = 300)

            # km.age_stat3 <- km.data |>
            #     mutate(age57_group = factor(ifelse(Age_dx > 57, "Age > 57", "Age <= 57"),levels = c("Age <= 57", "Age > 57"))) %>% 
            #     mutate(age57_stat3 = factor(paste(age57_group, STAT3_status, sep = " / "), levels = c("Age <= 57 / WT", "Age <= 57 / MT", "Age > 57 / WT", "Age > 57 / MT"))) |>
            #     dplyr::select(time_to_event, ever_tx_num, age57_stat3) |> na.omit()

            # fit.age_stat3 <- survfit(Surv(time_to_event, ever_tx_num) ~ age57_stat3, data = km.age_stat3)
    
            # p.km.age_stat3 <- ggsurvplot(
            #     fit.age_stat3, data = km.age_stat3,
            #     palette     = c("black", "#D55E00", "grey60", "#E69F00"),
            #     legend.labs = levels(km.age_stat3$age57_stat3),
            #     xlab        = "Time from diagnosis (months)",
            #     ylab        = "Treatment-free probability",
            #     pval        = TRUE, pval.size = 3,
            #     conf.int = TRUE, conf.int.alpha = 0.15,
            #     risk.table  = FALSE,
            #     ggtheme     = my_theme
            #     )$plot & theme_classic() + 
            #     theme(legend.position = c(0.8,0.7), 
            #         legend.title = element_blank(), 
            #         legend.text = element_text(size = 7),
            #         legend.key.height = unit(5, "mm"),
            #         legend.key.width  = unit(5, "mm"),

            #         legend.spacing.y = unit(0.5, "mm"),
            #         legend.box.spacing = unit(1, "mm"),
            #         legend.margin = margin(t = 0, r = 0, b = 0, l = 1))
        
            # ggsave(file.path(fig.path, "FigureS7.png"), p.km.age_stat3, width = 8, height = 8, units = "cm", dpi = 300)

        #     print(fit.age_stat3)


## overall survival analysis 
    
        os.dat <- cli_qced |>
            filter(time_to_treatment_qc_pass == TRUE, !is.na(Age_dx), !is.na(Sex)) |>
            mutate(
                dx_date_clean      = as.Date(dx_date_clean, format = "%m/%d/%y"),
                last_fu_raw_clean  = as.Date(last_fu_raw_clean, format = "%m/%d/%y"),
                deceased_raw_clean = as.Date(deceased_raw_clean, format = "%m/%d/%y"),
                os_event = case_when(
                Deceased_indicator == 2L             ~ 1L,
                Deceased_indicator %in% c(1L, 3L)   ~ 0L,
                TRUE                                 ~ NA_integer_
                ),
                os_end = case_when(
                Deceased_indicator == 2L & !is.na(deceased_raw_clean) ~ deceased_raw_clean,
                !is.na(last_fu_raw_clean)                             ~ last_fu_raw_clean,
                TRUE                                                  ~ NA_Date_
                ),
                os_time = as.numeric(os_end - dx_date_clean) / 30.44
            ) |> 
            filter(!is.na(os_event), !is.na(os_time), os_time > 0)

        os.cox <- coxph(Surv(os_time, os_event) ~ Age_dx + Sex, data = os.dat)
        print(summary(os.cox))
        print(cox.zph(os.cox))

## ── Figure 2; OS by age group  ────────────────────────────────────────────
    os.dat <- os.dat |>
    mutate(age60_group = factor(ifelse(Age_dx > 60, "Age > 60", "Age <= 60"),
                                levels = c("Age <= 60", "Age > 60")))

    km.os_age  <- os.dat |> dplyr::select(os_time, os_event, age60_group) |> na.omit()
    fit.os_age <- survfit(Surv(os_time, os_event) ~ age60_group, data = km.os_age)
    
    p.km.os_age <- ggsurvplot(
        fit.os_age, data = km.os_age,
        palette     = c("black", "grey60"),
        legend.title = "Age group",
        legend.labs = levels(km.os_age$age60_group),
        xlab        = "Time from diagnosis (months)",
        ylab        = "Overall survival probability",
        pval        = TRUE, pval.size = 3,
        conf.int = TRUE, conf.int.alpha = 0.15,
        risk.table  = FALSE,
        ggtheme     = my_theme
        )$plot + 
        theme_classic() + 
        theme(legend.position = c(0.2,.4), 
            legend.title = element_blank(), 
            legend.text = element_text(size = 7),
            legend.key.height = unit(5, "mm"),
            legend.key.width  = unit(5, "mm"),

            legend.spacing.y = unit(0.5, "mm"),
            legend.box.spacing = unit(1, "mm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 1))

    ggsave(file.path(fig.path, "Figure2B.png"), p.km.os_age, width = 8, height = 8, units = "cm", dpi = 300)

    # ggsave("figure2.tiff",
    #     (p.km.stat3 + labs(tag = "A")) | (p.km.os_age + labs(tag = "B")),
    #     width = two_column_width, height = 8, units = "cm", dpi = 300, compression = "lzw")
    

        os.cut.dat <- os.dat |>
        dplyr::select(os_time, os_event, Age_dx) |>
        filter(
            os_time > 0,
            os_event %in% c(0, 1)
        )

        ## Find the optimal age cutpoint
        ## minprop = 0.20: each group must contain at least 20% of patients
        cut.age <- surv_cutpoint(
            data      = os.cut.dat,
            time      = "os_time",
            event     = "os_event",
            variables = "Age_dx",
            minprop   = 0.20
        )

        print(cut.age)
        summary(cut.age)
