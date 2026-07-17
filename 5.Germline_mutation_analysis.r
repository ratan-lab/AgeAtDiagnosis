## Germline analysis for T-LGL leukemia
##
## Produces:
##   figure5.tiff   — CCR2 3'-UTR variant: age at diagnosis by carrier status (A)
##                    + Sanger trace placeholder (B)
##   figureS13.tiff — CCR2 3'-UTR RNA secondary structure arc diagram
##
## Individual-level germline genotypes are not distributed publicly (patient
## privacy). Two pre-computed summary files are provided instead:
##
##   germline_patient_summary.tsv.gz  — anonymised per-patient summary:
##     patient_id (anonymous), age, ccr2_carrier (0/1), lp_burden (total LP/P
##     variant count). Patient IDs are randomly shuffled and cannot be mapped
##     back to clinical identifiers.
##
##   germline_variant_summary.tsv.gz  — per-variant annotation table:
##     variant_key, gene, func, ACMG_prediction, n_carriers, pct_carriers.
##     Variant-level Firth association tests (logistf, firth=TRUE) were run on
##     the individual-level data and are not re-computable from these summaries;
##     results are reported in Supplementary Table S5.
##
## Variant-level Firth analysis strategy (documented here for methods
## transparency):
##   - LP/P variants filtered to LGL-biology gene universe (driver genes +
##     Hallmark IL2/STAT5, inflammatory, apoptosis, PI3K/AKT, IFN-γ;
##     KEGG JAK-STAT, T-cell, NK-cell, apoptosis; manually curated T-LGL genes)
##   - Variants with >= 3 carriers and carrier rate <= 20% tested
##   - Firth penalised logistic: carrier ~ Age (logistf, firth = TRUE)
##   - BH correction across all tested variants

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(Biostrings)
  library(GenomicRanges)
  library(tidyverse)
  library(logistf)
  library(patchwork)
  library(httr)
  library(jsonlite)
})

source("funs.R")

MIN_ALT <- 3

## ── Load summary data ─────────────────────────────────────────────────────────
germ_patients <- read_tsv(germline.patients.path, show_col_types = FALSE)
germ_variants <- read_tsv(germline.variants.path, show_col_types = FALSE)

cat(sprintf("Germline patient summary: %d patients\n", nrow(germ_patients)))
cat(sprintf("Germline variant summary: %d unique LP/P variants\n", nrow(germ_variants)))
cat(sprintf("Variants with >= %d carriers: %d\n",
            MIN_ALT, sum(germ_variants$n_carriers >= MIN_ALT, na.rm = TRUE)))

## ── Germline burden vs age ────────────────────────────────────────────────────
sp_burden <- cor.test(germ_patients$age, germ_patients$lp_burden,
                      method = "spearman", exact = FALSE)
cat(sprintf("\nGermline LP/P burden vs Age: rho = %.3f, p = %.4f\n",
            sp_burden$estimate, sp_burden$p.value))

## ── CCR2 3'-UTR variant focused analysis ─────────────────────────────────────
## CCR2 chr3:46400346 C>T (rs536776017) — classified LP by ACMG-tapes
## (Prob_Path = 0.32); validated by Sanger sequencing.

ccr2_df    <- germ_patients |> mutate(CCR2_carrier = as.logical(ccr2_carrier))
n_ccr2     <- sum(ccr2_df$ccr2_carrier, na.rm = TRUE)
n_total    <- nrow(ccr2_df)

cat(sprintf("\nCCR2 UTR3 carriers: %d / %d (%.1f%%)\n",
            n_ccr2, n_total, 100 * n_ccr2 / n_total))

age_car    <- na.omit(ccr2_df$age[ccr2_df$CCR2_carrier])
age_noncar <- na.omit(ccr2_df$age[!ccr2_df$CCR2_carrier])
wt_ccr2    <- wilcox.test(age_car, age_noncar)

cat(sprintf("Age: carriers median = %.1f (n=%d), non-carriers median = %.1f (n=%d)\n",
            median(age_car), length(age_car), median(age_noncar), length(age_noncar)))
cat(sprintf("Wilcoxon p = %.4f\n", wt_ccr2$p.value))

ok_ccr2   <- !is.na(ccr2_df$age)
m_ccr2    <- logistf(ccr2_df$CCR2_carrier[ok_ccr2] ~ ccr2_df$age[ok_ccr2], firth = TRUE)
age_pname <- grep("age", names(m_ccr2$prob), value = TRUE)[1]
firth_p   <- m_ccr2$prob[age_pname]
cat(sprintf("Firth logistic: age coef = %.4f, p = %.4f (n=%d carriers; interpret with caution)\n",
            m_ccr2$coefficients[age_pname], firth_p, n_ccr2))

## ── Variant-level summary (Firth tests pre-computed; see header note) ─────────
cat("\nTop LP/P variants by carrier frequency in LGL cohort:\n")
print(germ_variants |>
        filter(n_carriers >= MIN_ALT) |>
        arrange(desc(n_carriers)) |>
        dplyr::select(variant_key, Gene, Func, ACMG_prediction, n_carriers, pct_carriers) |>
        head(20), n = 20)

## ── LGL-biology gene universe (documented for methods transparency) ────────────
driver_genes <- c("STAT3","STAT5B","TET2","DNMT3A","KMT2D","KDM6A","TNFAIP3",
                  "PIK3R1","ABCC9","FAS","SULF1","SLC6A15","DDX59","ARHGAP25","PCDHA11")
extra_genes  <- c(
  "JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT4","STAT5A","STAT6",
  "SH2B3","SOCS1","SOCS3","IL2","IL2RA","IL2RB","IL2RG","IL15","IL15RA",
  "PRF1","GZMB","GZMM","GZMH","GZMK","GNLY","NKG7",
  "KLRD1","KLRC1","KLRC2","KLRB1","KLRK1","NCR1","NCR3","CD226",
  "BCL2","BCL2L1","MCL1","BCL2L11","BID","CASP3","CASP8","CFLAR",
  "PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT",
  "CCR1","CCR2","CCR3","CCR4","CCR5","CCR7","CXCR3","CXCR4","CXCR6",
  "TET1","TET3","DNMT1","DNMT3B","KMT2A","EZH2","ASXL1","ASXL2"
)
## Full gene universe additionally includes Hallmark and KEGG gene sets
## (retrieved via msigdbr in the primary analysis); see manuscript Methods.

cat(sprintf("\nLGL driver + curated genes: %d\n",
            length(unique(c(driver_genes, extra_genes)))))

## ── RNAfold analysis — CCR2 3'-UTR secondary structure ───────────────────────
FLANK   <- 50L
CCR2_CHR  <- "chr3"
CCR2_POS  <- 46400346L
CCR2_REF  <- "C"
CCR2_ALT  <- "T"
CCR2_RSID <- "rs536776017"

genome     <- BSgenome.Hsapiens.UCSC.hg19
context_gr <- GRanges(CCR2_CHR, IRanges(CCR2_POS - FLANK, CCR2_POS + FLANK))
seq_genomic_ref <- as.character(getSeq(genome, context_gr))
seq_genomic_alt <- seq_genomic_ref
substr(seq_genomic_alt, FLANK + 1L, FLANK + 1L) <- CCR2_ALT

rna_ref <- gsub("T", "U", seq_genomic_ref)
rna_alt <- gsub("T", "U", seq_genomic_alt)

RNAFOLD_BIN <- Sys.which("RNAfold")
if (nchar(RNAFOLD_BIN) == 0)
  stop("RNAfold not found in PATH. Install via: conda install -c bioconda viennarna")

run_rnafold <- function(seq, label = "") {
  tmp_in  <- tempfile(fileext = ".fa")
  tmp_out <- tempfile(fileext = ".txt")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)
  writeLines(c(paste0(">", label), seq), tmp_in)
  rc <- system2(RNAFOLD_BIN,
                args   = c("--noPS", "--dangles=2", "--temp=37"),
                stdin  = tmp_in, stdout = tmp_out, stderr = FALSE)
  if (rc != 0) stop("RNAfold exited with code ", rc)
  lines       <- readLines(tmp_out)
  struct_line <- lines[grepl("^[.()+&]+\\s+\\(", lines)][1]
  if (is.na(struct_line))
    stop("RNAfold produced no parseable structure for label '", label, "'")
  structure <- strsplit(trimws(struct_line), "\\s+")[[1]][1]
  mfe       <- as.numeric(gsub("[() ]", "", sub("^\\S+\\s+", "", trimws(struct_line))))
  list(label = label, sequence = seq, structure = structure, mfe = mfe)
}

fold_ref <- run_rnafold(rna_ref, "Reference_C")
fold_alt <- run_rnafold(rna_alt, "Alternate_T")
ddG      <- fold_alt$mfe - fold_ref$mfe

cat(sprintf("Reference MFE : %.2f kcal/mol\n", fold_ref$mfe))
cat(sprintf("Alternate MFE : %.2f kcal/mol\n", fold_alt$mfe))
cat(sprintf("ΔΔMFE (alt-ref): %.2f kcal/mol\n", ddG))

n_nt      <- nchar(fold_ref$structure)
ref_chars <- strsplit(fold_ref$structure, "")[[1]]
alt_chars <- strsplit(fold_alt$structure, "")[[1]]
cat(sprintf("Positions with changed pairing: %d / %d\n",
            sum(ref_chars != alt_chars), n_nt))

## ── RegulomeDB ────────────────────────────────────────────────────────────────
regulome_resp <- tryCatch(
  GET("https://regulomedb.org/regulome-search/",
      query   = list(regions = "chr3:46400345-46400346",
                     genome  = "GRCh37", format = "json"),
      timeout(30)),
  error = function(e) { message("RegulomeDB request failed: ", e$message); NULL })

if (!is.null(regulome_resp) && !http_error(regulome_resp)) {
  reg_json <- tryCatch(fromJSON(content(regulome_resp, "text", encoding = "UTF-8"),
                                flatten = TRUE), error = function(e) NULL)
  if (!is.null(reg_json)) {
    score_val <- reg_json$score %||% reg_json$regulome_score %||%
                 tryCatch(reg_json$`@graph`[[1]]$score, error = function(e) NULL)
    if (!is.null(score_val)) cat(sprintf("RegulomeDB score: %s\n", score_val))
    features <- reg_json$features
    if (!is.null(features) && length(features) > 0) {
      cat("RegulomeDB features:\n"); print(as.data.frame(features))
    }
  }
} else {
  cat("RegulomeDB query unavailable\n")
}

################################################################################
## Manuscript figures
################################################################################

## Figure 5 — CCR2 variant: age by carrier status (A) + Sanger placeholder (B)
n_lab <- sprintf("Carrier\n(n=%d)", n_ccr2)

panel_a <- ccr2_df |>
  mutate(CCR2_status = ifelse(CCR2_carrier, n_lab, "Non-carrier"),
         CCR2_status = factor(CCR2_status, levels = c("Non-carrier", n_lab))) |>
  ggplot(aes(CCR2_status, age, colour = CCR2_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, colour = "grey40") +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.2) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"), guide = "none") +
  annotate("text", x = 1.5, y = max(ccr2_df$age, na.rm = TRUE) + 3,
           label = sprintf("Wilcoxon p = %.3f\nFirth p = %.3f", wt_ccr2$p.value, firth_p),
           size = 2.5, hjust = 0.5) +
  labs(x = NULL, y = "Age at diagnosis (years)") +
  my_theme

panel_b <- ggplot() +
  annotate("rect", xmin = 0.05, xmax = 0.95, ymin = 0.1, ymax = 0.9,
           fill = "grey93", colour = "grey70", linewidth = 0.4) +
  annotate("text", x = 0.5, y = 0.65,
           label = "Sanger sequencing trace",
           size = 2.8, colour = "grey45", fontface = "italic") +
  annotate("text", x = 0.5, y = 0.50,
           label = "CCR2 3'-UTR  chr3:46400346 C>T",
           size = 2.3, colour = "grey45") +
  annotate("text", x = 0.5, y = 0.35,
           label = "(chromatogram placeholder)",
           size = 2.0, colour = "grey60") +
  theme_void() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave("figure5.tiff",
       panel_a + panel_b +
         plot_annotation(tag_levels = "A") &
         theme(plot.tag = element_text(size = 9, face = "bold")),
       width = 14, height = 8, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figure5.tiff\n")

## Figure S13 — CCR2 3'-UTR RNA secondary structure arc diagram
parse_pairs <- function(structure) {
  chars <- strsplit(structure, "")[[1]]
  stack <- integer(0)
  pairs <- data.frame(i = integer(), j = integer())
  for (k in seq_along(chars)) {
    if (chars[k] == "(") {
      stack <- c(stack, k)
    } else if (chars[k] == ")") {
      if (length(stack) == 0) next
      i     <- stack[length(stack)]
      stack <- stack[-length(stack)]
      pairs <- rbind(pairs, data.frame(i = i, j = k))
    }
  }
  pairs
}

make_arcs <- function(pairs, allele, n_nt, variant_pos) {
  if (nrow(pairs) == 0) return(tibble())
  pmap_dfr(pairs, function(i, j) {
    centre <- (i + j) / 2
    radius <- (j - i) / 2
    theta  <- seq(0, pi, length.out = 40)
    tibble(
      x             = centre + radius * cos(theta),
      y             = radius * sin(theta),
      arc_id        = paste(allele, i, j, sep = "_"),
      allele        = allele,
      spans_variant = (i <= variant_pos & j >= variant_pos)
    )
  })
}

pairs_ref <- parse_pairs(fold_ref$structure)
pairs_alt <- parse_pairs(fold_alt$structure)
arc_all   <- bind_rows(
  make_arcs(pairs_ref, "Reference (C)", n_nt, FLANK + 1L),
  make_arcs(pairs_alt, "Alternate (T)", n_nt, FLANK + 1L)
) |> mutate(allele = factor(allele, levels = c("Reference (C)", "Alternate (T)")))

seq_bar <- bind_rows(
  tibble(pos = seq_len(n_nt), base = strsplit(rna_ref, "")[[1]],
         is_var = pos == (FLANK + 1L), allele = "Reference (C)"),
  tibble(pos = seq_len(n_nt), base = strsplit(rna_alt, "")[[1]],
         is_var = pos == (FLANK + 1L), allele = "Alternate (T)")
) |> mutate(allele = factor(allele, levels = c("Reference (C)", "Alternate (T)")))

p_arc <- ggplot() +
  geom_path(data = arc_all,
            aes(x, y, group = arc_id,
                colour = spans_variant, linewidth = spans_variant), alpha = 0.75) +
  geom_point(data = seq_bar,
             aes(x = pos, y = 0, fill = is_var),
             shape = 21, size = 1.2, colour = "grey30", stroke = 0.2) +
  scale_colour_manual(values = c("FALSE" = "#56B4E9", "TRUE" = "#E69F00"),
                      name = "Spans variant") +
  scale_linewidth_manual(values = c("FALSE" = 0.4, "TRUE" = 1.0),
                         name = "Spans variant") +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#E69F00"),
                    name = "Variant site") +
  facet_wrap(~allele, ncol = 1) +
  annotate("text", x = FLANK + 1L, y = -2, label = "▲", size = 3, colour = "#E69F00") +
  expand_limits(y = -2.5) +
  labs(x = "Position in 101 nt window (nt 245 of CCR2 3'UTR at centre)",
       y = "Arc height (proportional to stem span)",
       title = sprintf("CCR2 3'-UTR secondary structure: %s", CCR2_RSID),
       subtitle = sprintf("Ref MFE = %.2f  |  Alt MFE = %.2f  |  ΔΔMFE = %.2f kcal/mol",
                          fold_ref$mfe, fold_alt$mfe, ddG)) +
  my_theme +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold", size = 9))

ggsave("figureS13.tiff", p_arc,
       width = 17, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS13.tiff\n")
