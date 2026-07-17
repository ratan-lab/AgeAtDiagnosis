## Cohort summary statistics for Table 1
##
## Produces: console output summarised in Table 1 of the manuscript

suppressPackageStartupMessages(library(tidyverse))
source("funs.R")

cli <- read_csv(clinical.path, show_col_types = FALSE)

cat(sprintf("N = %d\n", nrow(cli)))

## ── Sex ──────────────────────────────────────────────────────────────────────
sex_tab <- cli |> count(Sex)
f <- filter(sex_tab, Sex == "F")$n
cat(sprintf("Female: %d (%.1f%%)\n", f, 100 * f / nrow(cli)))

## ── Age at diagnosis ─────────────────────────────────────────────────────────
ages <- cli$Age_dx
cat(sprintf("Age at dx: median %.0f, IQR (%.0f, %.0f), range %.0f–%.0f, n = %d\n",
            median(ages, na.rm = TRUE),
            quantile(ages, 0.25, na.rm = TRUE),
            quantile(ages, 0.75, na.rm = TRUE),
            min(ages, na.rm = TRUE),
            max(ages, na.rm = TRUE),
            sum(!is.na(ages))))

## ── RA ───────────────────────────────────────────────────────────────────────
ra_known <- cli |> filter(RA %in% c("Y", "N"))
ra_y <- sum(ra_known$RA == "Y")
cat(sprintf("RA+: %d/%d (%.1f%%) [excludes %d with unknown RA status]\n",
            ra_y, nrow(ra_known), 100 * ra_y / nrow(ra_known),
            nrow(cli) - nrow(ra_known)))

## ── STAT3 ────────────────────────────────────────────────────────────────────
stat3_known <- cli |> filter(STAT3_status %in% c("MT", "WT"))
stat3_mt <- sum(stat3_known$STAT3_status == "MT")
cat(sprintf("STAT3 MT: %d/%d (%.1f%%) [excludes %d untested]\n",
            stat3_mt, nrow(stat3_known), 100 * stat3_mt / nrow(stat3_known),
            nrow(cli) - nrow(stat3_known)))

## ── CBC at diagnosis ─────────────────────────────────────────────────────────
cbc_vars <- c(WBC = "WBC", HGB = "HGB", ANC = "ANC", ALC = "ALC", Platelets = "Platelets")

for (label in names(cbc_vars)) {
  col <- cbc_vars[label]
  vals <- as.numeric(cli[[col]])
  vals <- vals[!is.na(vals)]
  cat(sprintf("%-10s  median %.2f, IQR (%.2f, %.2f), n = %d\n",
              label,
              median(vals),
              quantile(vals, 0.25),
              quantile(vals, 0.75),
              length(vals)))
}

## ── Treatment ────────────────────────────────────────────────────────────────
tx_known <- cli |> filter(ever_tx %in% c("Y", "N"))
tx_y <- sum(tx_known$ever_tx == "Y")
cat(sprintf("\nEver treated: %d/%d (%.1f%%)\n",
            tx_y, nrow(tx_known), 100 * tx_y / nrow(tx_known)))

# First-line agent breakdown (among treated)
tx1_counts <- cli |>
  filter(ever_tx == "Y", !is.na(Tx1), Tx1 != "NA") |>
  count(Tx1, sort = TRUE)
cat("First-line agent (treated patients):\n")
print(tx1_counts, n = 10)

# Response (among treated with evaluable response)
resp_eval <- cli |>
  filter(ever_tx == "Y") |>
  mutate(resp = toupper(trimws(Tx1_response))) |>
  filter(resp %in% c("CR", "PR", "NR"))
cr <- sum(resp_eval$resp == "CR")
pr <- sum(resp_eval$resp == "PR")
nr <- sum(resp_eval$resp == "NR")
cat(sprintf("Treatment response (n = %d evaluable): CR = %d, PR = %d, NR = %d\n",
            nrow(resp_eval), cr, pr, nr))
cat(sprintf("Overall response rate (CR+PR): %d/%d (%.1f%%)\n",
            cr + pr, nrow(resp_eval), 100 * (cr + pr) / nrow(resp_eval)))

## ── Follow-up / time-to-event (QC-passed subset) ─────────────────────────────
qc <- cli |> filter(time_to_treatment_qc_pass == TRUE)
cat(sprintf("\nQC-passed for time-to-treatment analysis: %d\n", nrow(qc)))
cat(sprintf("  Treated: %d  |  Not treated / censored: %d\n",
            sum(qc$ever_tx == "Y", na.rm = TRUE),
            sum(qc$ever_tx == "N", na.rm = TRUE)))
tte <- qc$time_to_event[!is.na(qc$time_to_event)]
cat(sprintf("  Median time to event: %.1f months, IQR (%.1f, %.1f)\n",
            median(tte), quantile(tte, 0.25), quantile(tte, 0.75)))

## ── Genomic subsets ──────────────────────────────────────────────────────────
n_wes <- sum(!is.na(cli$patient_id) & cli$patient_id != "")
cat(sprintf("\nWES + germline data available: %d\n", n_wes))
