## Time-to-treatment and survival analysis for T-LGL leukemia
##
## Produces:
##   figure2.tiff  — STAT3 KM treatment-free survival (Fig 2)
##   figureS5.tiff — Age × STAT3 four-group KM (Fig S5)
##   figureS6.tiff — RA KM (Fig S6)

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(splines)
  library(cmprsk)
  library(logistf)
  library(patchwork)
})

source("funs.R")

## ── Load data ─────────────────────────────────────────────────────────────────
cli_qced <- read_csv(clinical.path, show_col_types = FALSE) |>
  mutate(
    Age_dx                    = as.numeric(Age_dx),
    Sex                       = factor(Sex, levels = c("F", "M")),
    STAT3_status              = factor(STAT3_status, levels = c("WT", "MT")),
    time_to_treatment_qc_pass = as.logical(time_to_treatment_qc_pass),
    death_used_for_censor     = as.logical(death_used_for_censor)
  )

cat("Total T-LGL N =", nrow(cli_qced), "\n")
cat("QC passed =", sum(cli_qced$time_to_treatment_qc_pass, na.rm = TRUE), "\n")
cat("Treated events =", sum(cli_qced$time_to_treatment_qc_pass & cli_qced$ever_tx == "Y",
                            na.rm = TRUE), "\n")
cat("Median time to event (months) =",
    median(cli_qced$time_to_event[cli_qced$time_to_treatment_qc_pass], na.rm = TRUE), "\n")

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
                Age_dx, Sex, STAT3_status, RA_status) |>
  filter(!is.na(ever_tx_num), !is.na(Age_dx), !is.na(Sex),
         !is.na(STAT3_status), !is.na(RA_status))

cat("Analytic N =", nrow(cox.master), "| Events (treated) =",
    sum(cox.master$ever_tx_num == 1), "\n")

################################################################################
## COX MODELS
################################################################################

## ── Main model: Age + Sex ─────────────────────────────────────────────────────
main.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex, data = cox.master)
print(summary(main.cox))
print(cox.zph(main.cox))

## ── Non-linear age test (LRT: linear vs natural cubic spline) ────────────────
cox.age.spline <- coxph(Surv(time_to_event, ever_tx_num) ~ ns(Age_dx, df = 3) + Sex,
                        data = cox.master)
cat("\nLRT: linear vs cubic spline for Age_dx\n")
print(anova(main.cox, cox.age.spline))

## ── STAT3 model ──────────────────────────────────────────────────────────────
stat3.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + STAT3_status,
                   data = cox.master)
print(summary(stat3.cox))
print(cox.zph(stat3.cox))

## ── RA model ─────────────────────────────────────────────────────────────────
ra.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + RA_status,
                data = cox.master)
print(summary(ra.cox))

## ── Fully adjusted model ─────────────────────────────────────────────────────
full.cox <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx + Sex + STAT3_status + RA_status,
                  data = cox.master)
print(summary(full.cox))

## ── Interaction models ────────────────────────────────────────────────────────
stat3.int <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx * STAT3_status + Sex,
                   data = cox.master)
cat("\nLRT: Age × STAT3 interaction\n")
print(anova(stat3.cox, stat3.int, test = "LRT"))

ra.int <- coxph(Surv(time_to_event, ever_tx_num) ~ Age_dx * RA_status + Sex,
                data = cox.master)
cat("\nLRT: Age × RA interaction\n")
print(anova(ra.cox, ra.int, test = "LRT"))

################################################################################
## COMPETING RISKS (Fine-Gray)
################################################################################
## cr_status: 0 = censored, 1 = treatment, 2 = death before treatment

cat("\nCompeting-risk event counts:\n")
print(table(cox.master$cr_status, dnn = "0=censored/1=treatment/2=death"))

n_cr <- nrow(cox.master)
n_tx <- sum(cox.master$cr_status == 1)

cr.covs.main  <- model.matrix(~ Age_dx + Sex, data = cox.master)[, -1, drop = FALSE]
fg.main       <- crr(ftime = cox.master$time_to_event, fstatus = cox.master$cr_status,
                     cov1 = cr.covs.main, failcode = 1, cencode = 0)
cat("\nFine-Gray: Age_dx + Sex\n"); print(summary(fg.main))

cr.covs.stat3 <- model.matrix(~ Age_dx + Sex + STAT3_status, data = cox.master)[, -1, drop = FALSE]
fg.stat3      <- crr(ftime = cox.master$time_to_event, fstatus = cox.master$cr_status,
                     cov1 = cr.covs.stat3, failcode = 1, cencode = 0)
cat("\nFine-Gray: Age_dx + Sex + STAT3_status\n"); print(summary(fg.stat3))

cr.covs.ra <- model.matrix(~ Age_dx + Sex + RA_status, data = cox.master)[, -1, drop = FALSE]
fg.ra      <- crr(ftime = cox.master$time_to_event, fstatus = cox.master$cr_status,
                  cov1 = cr.covs.ra, failcode = 1, cencode = 0)
cat("\nFine-Gray: Age_dx + Sex + RA_status\n"); print(summary(fg.ra))

################################################################################
## KAPLAN-MEIER CURVES
################################################################################

km.data <- cox.master |>
  mutate(age60_group = factor(ifelse(Age_dx > 60, "Age > 60", "Age <= 60"),
                               levels = c("Age <= 60", "Age > 60")))

## ── STAT3 (→ Figure 2) ───────────────────────────────────────────────────────
km.stat3  <- km.data |> dplyr::select(time_to_event, ever_tx_num, STAT3_status) |> na.omit()
fit.stat3 <- survfit(Surv(time_to_event, ever_tx_num) ~ STAT3_status, data = km.stat3)
print(fit.stat3)
print(survdiff(Surv(time_to_event, ever_tx_num) ~ STAT3_status, data = km.stat3))

## ── RA (→ Figure S6) ─────────────────────────────────────────────────────────
km.ra  <- km.data |> dplyr::select(time_to_event, ever_tx_num, RA_status) |> na.omit()
fit.ra <- survfit(Surv(time_to_event, ever_tx_num) ~ RA_status, data = km.ra)
print(fit.ra)
print(survdiff(Surv(time_to_event, ever_tx_num) ~ RA_status, data = km.ra))

## ── Age × STAT3 four-group (→ Figure S5) ─────────────────────────────────────
km.age_stat3 <- km.data |>
  mutate(age60_stat3 = factor(
    paste(age60_group, STAT3_status, sep = " / "),
    levels = c("Age <= 60 / WT", "Age <= 60 / MT", "Age > 60 / WT", "Age > 60 / MT")
  )) |>
  dplyr::select(time_to_event, ever_tx_num, age60_stat3) |>
  na.omit()

fit.age_stat3 <- survfit(Surv(time_to_event, ever_tx_num) ~ age60_stat3, data = km.age_stat3)
print(fit.age_stat3)

################################################################################
## OVERALL SURVIVAL
################################################################################

os.dat <- cli_qced |>
  filter(time_to_treatment_qc_pass == TRUE, !is.na(Age_dx), !is.na(Sex)) |>
  mutate(
    dx_date_clean      = as.Date(dx_date_clean),
    last_fu_raw_clean  = as.Date(last_fu_raw_clean),
    deceased_raw_clean = as.Date(deceased_raw_clean),
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

cat(sprintf("\nOS analysis: N = %d, deaths = %d\n", nrow(os.dat), sum(os.dat$os_event)))

os.cox <- coxph(Surv(os_time, os_event) ~ Age_dx + Sex, data = os.dat)
print(summary(os.cox))
print(cox.zph(os.cox))

## ── OS by age group (→ Figure 2B) ────────────────────────────────────────────
os.dat <- os.dat |>
  mutate(age60_group = factor(ifelse(Age_dx > 60, "Age > 60", "Age <= 60"),
                               levels = c("Age <= 60", "Age > 60")))

km.os_age  <- os.dat |> dplyr::select(os_time, os_event, age60_group) |> na.omit()
fit.os_age <- survfit(Surv(os_time, os_event) ~ age60_group, data = km.os_age)
print(fit.os_age)
print(survdiff(Surv(os_time, os_event) ~ age60_group, data = km.os_age))

################################################################################
## TREATMENT RESPONSE
################################################################################

resp.dat <- cli_qced |>
  filter(time_to_treatment_qc_pass == TRUE, !is.na(Age_dx),
         ever_tx == "Y", !is.na(Tx1_response)) |>
  mutate(
    resp_clean = toupper(trimws(Tx1_response)),
    responder  = case_when(
      resp_clean %in% c("CR", "PR") ~ 1L,
      resp_clean == "NR"            ~ 0L,
      TRUE                          ~ NA_integer_
    )
  ) |>
  filter(!is.na(responder))

cat(sprintf("Treatment response: N = %d, R = %d, NR = %d\n",
            nrow(resp.dat), sum(resp.dat$responder), sum(resp.dat$responder == 0)))

resp.firth <- logistf(responder ~ Age_dx, data = resp.dat)
cat(sprintf("Treatment response Firth p (Age_dx) = %.4f\n", resp.firth$prob["Age_dx"]))

################################################################################
## MANUSCRIPT FIGURES
################################################################################

## figure2.tiff  — panel A: STAT3 KM | panel B: OS by age group
## figureS5.tiff — Age × STAT3 four-group KM (Fig S5)
## figureS6.tiff — RA KM (Fig S6)

p.km.stat3 <- ggsurvplot(
  fit.stat3, data = km.stat3,
  palette     = c("black", "grey60"),
  legend.labs = levels(km.stat3$STAT3_status),
  xlab        = "Time from diagnosis (months)",
  ylab        = "Treatment-free probability",
  pval        = TRUE, pval.size = 3,
  conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table  = FALSE,
  ggtheme     = my_theme
)$plot

p.km.os_age <- ggsurvplot(
  fit.os_age, data = km.os_age,
  palette     = c("black", "grey60"),
  legend.labs = levels(km.os_age$age60_group),
  xlab        = "Time from diagnosis (months)",
  ylab        = "Overall survival probability",
  pval        = TRUE, pval.size = 3,
  conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table  = FALSE,
  ggtheme     = my_theme
)$plot

ggsave("figure2.tiff",
       (p.km.stat3 + labs(tag = "A")) | (p.km.os_age + labs(tag = "B")),
       width = two_column_width, height = 8, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figure2.tiff\n")

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
)$plot

ggsave("figureS5.tiff", p.km.age_stat3,
       width = onepointfive_column_width, height = 8, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS5.tiff\n")

p.km.ra <- ggsurvplot(
  fit.ra, data = km.ra,
  palette     = c("black", "grey60"),
  legend.labs = levels(km.ra$RA_status),
  xlab        = "Time from diagnosis (months)",
  ylab        = "Treatment-free probability",
  pval        = TRUE, pval.size = 3,
  conf.int = TRUE, conf.int.alpha = 0.15,
  risk.table  = FALSE,
  ggtheme     = my_theme
)$plot

ggsave("figureS6.tiff", p.km.ra,
       width = onepointfive_column_width, height = 8, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS6.tiff\n")
