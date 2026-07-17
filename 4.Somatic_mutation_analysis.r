## This script is for somatic mutation analysis and creating the plots

### Set Up Packages and paths
suppressPackageStartupMessages({
library(ComplexHeatmap)
library(tidyverse)
library(data.table)
library(logistf)
})

source("funs.R")

## ── Clinical data ──────────────────────────────────────────────────────────
cli <- as.data.frame(fread(clinical.path)) |>
    rename(Age = Age_dx) |>
    filter(!is.na(patient_id) & patient_id != "") |>
    mutate(RegID.2 = paste0(patient_id, "_tumor"))

## ── Somatic mutation data ───────────────────────────────────────────────────
maf <- read_tsv(somatic.path, show_col_types = FALSE)

# gene × sample count matrix
mut.mat <- maf |>
    count(Hugo_Symbol, Tumor_Sample_Barcode) |>
    pivot_wider(names_from = Tumor_Sample_Barcode, values_from = n,
                values_fill = 0L) |>
    column_to_rownames("Hugo_Symbol") |>
    as.matrix()

# Subset overlapped samples
common.ids   <- intersect(cli$RegID.2, colnames(mut.mat))
cli.filt     <- cli[match(common.ids, cli$RegID.2), ]
mut.mat.filt <- mut.mat[, match(common.ids, colnames(mut.mat))]

# capture TMB (raw counts) before binarisation
tmb_vec <- colSums(mut.mat.filt)

# prepare the gene matrix
mut.mat.filt[mut.mat.filt != 0] <- 1
mut.mat.filt <- mut.mat.filt[rowSums(mut.mat.filt) > 0, ]

## ── TMB vs age ──────────────────────────────────────────────────────────────

tmb_df <- tibble(
  ID           = colnames(mut.mat.filt),
  TMB          = tmb_vec,
  Age          = cli.filt$Age,
  STAT3_status = factor(cli.filt$STAT3_status, levels = c("WT", "MT"))
)

cor_tmb <- cor.test(tmb_df$Age, tmb_df$TMB, method = "spearman", exact = FALSE)
cat(sprintf("TMB vs Age: rho=%.3f, p=%.3g\n", cor_tmb$estimate, cor_tmb$p.value))

lm_adj    <- lm(TMB ~ Age + STAT3_status, data = filter(tmb_df, !is.na(STAT3_status)))
lm_coefs  <- coef(summary(lm_adj))
lm_age_p  <- lm_coefs["Age", "Pr(>|t|)"]
cat(sprintf("Adjusted lm TMB ~ Age + STAT3: Age coef=%.3f, p=%.4f\n",
            lm_coefs["Age", "Estimate"], lm_age_p))
print(summary(lm_adj))

## ── generalized linear model for each gene ──────────────────────────────────

# minimum number of mutated samples to include a gene/pathway in testing;
# pre-filtering before BH correction reduces the multiple-testing burden
MIN_MUT <- 5

mut_test      <- mut.mat.filt[rowSums(mut.mat.filt) >= MIN_MUT, ]
cat(sprintf("Testing %d genes (of %d) with >= %d mutations\n",
            nrow(mut_test), nrow(mut.mat.filt), MIN_MUT))

# Firth penalised logistic regression (guards against separation with rare events)
age = cli.filt$Age
df <- tibble()
for (g in 1:nrow(mut_test)) {
  y <- mut_test[g,]
  x <- age
  gene.name = rownames(mut_test)[g]
  model = logistf(y ~ 1 + x, firth = TRUE)
  p.logistf = model$prob[2]
  df <- bind_rows(df, tibble(gene = gene.name, pval.logistf = p.logistf,
                             num_samples = sum(y), coef = model$coefficients["x"]))
}

df.genes <- df |> mutate(padj = p.adjust(pval.logistf, method = "BH"))
print(df.genes |> arrange(padj))
sig.genes <- df.genes |> filter(padj < 0.05) |> pull(gene)

## ── Wilcoxon test: age distribution in mutated vs unmutated patients ────────
##    only in driven genes.

TLGL_DRIVER_GENES <- c("STAT3", "KMT2D", "TNFAIP3", "KDM6A", "ABCC9", "PIK3R1",
                       "TET2", "PCDHA11", "SLC6A15", "SULF1", "ARHGAP25", "DDX59",
                       "DNMT3A", "FAS", "STAT5B")
driver.present  <- intersect(TLGL_DRIVER_GENES, rownames(mut.mat.filt))

wilcox_res <- lapply(driver.present, function(g) {
  mut_status <- mut.mat.filt[g, ]
  age_mut    <- cli.filt$Age[mut_status == 1]
  age_unmut  <- cli.filt$Age[mut_status == 0]
  if (length(age_mut) < 2 || length(age_unmut) < 2) return(NULL)
  wt <- wilcox.test(age_mut, age_unmut)
  tibble(gene        = g,
         n_mutated   = length(age_mut),
         n_unmutated = length(age_unmut),
         median_mut  = median(age_mut),
         median_unmut= median(age_unmut),
         median_diff = median(age_mut) - median(age_unmut),
         pval        = wt$p.value)
}) |> bind_rows() |>
  mutate(padj = p.adjust(pval, method = "BH"))

cat("\nWilcoxon test — age (mutated vs unmutated) for driver genes:\n")
print(wilcox_res |> arrange(pval), n = 20)

## figure4.tiff — TMB vs age, colored by STAT3 status
cor_label <- sprintf("Spearman rho = %.2f, p = %.3g\nAge adj. for STAT3: p = %.3g",
                     cor_tmb$estimate, cor_tmb$p.value, lm_age_p)

figure4 <- ggplot(filter(tmb_df, !is.na(STAT3_status)),
                  aes(Age, TMB, color = STAT3_status, fill = STAT3_status)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  scale_color_manual(values = c(WT = "black", MT = "#D55E00"), name = "STAT3") +
  scale_fill_manual(values  = c(WT = "black", MT = "#D55E00"), guide = "none") +
  annotate("text", x = Inf, y = Inf, label = cor_label,
           hjust = 1.1, vjust = 1.5, size = 3) +
  labs(x = "Age at diagnosis (years)",
       y = "Tumor mutation burden") +
  my_theme

ggsave("figure4.tiff", figure4,
       width = onepointfive_column_width, height = 7, units = "cm", dpi = 300,
       compression = "lzw")

## figureS11.tiff — Driver gene OncoPrint sorted by age (14 x 6 cm)
##   restricted to driver genes with >= 5 mutated patients

s11.genes <- intersect(TLGL_DRIVER_GENES, names(which(rowSums(mut.mat.filt) >= 5)))
s11.mut   <- mut.mat.filt[s11.genes, ]

s11.mut <- matrix(
  ifelse(s11.mut == 1, "Mut", " "),
  nrow = nrow(s11.mut), ncol = ncol(s11.mut),
  dimnames = dimnames(s11.mut)
)

mut_colours <- c("Mut" = "#008000")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "azure2", col = NA))
  },
  Mut = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = mut_colours["Mut"], col = NA))
  }
)

pat.anno.dat <- data.frame(Age = cli.filt$Age, ID = cli.filt$RegID.2)
pat.anno.dat <- pat.anno.dat[order(pat.anno.dat$Age), ]
pat.anno <- pat.anno.dat$Age
names(pat.anno) <- pat.anno.dat$ID

col_order <- match(names(pat.anno), colnames(s11.mut))
stopifnot(!anyNA(col_order))
s11.mut.order <- s11.mut[, col_order]

age_mid  <- mean(pat.anno)
age_half <- max(abs(pat.anno - age_mid))

tiff("figureS11.tiff",
     width = onepointfive_column_width, height = 3, units = "cm", res = 300,
     compression = "lzw")
oncoPrint(s11.mut.order,
          alter_fun = alter_fun, col = mut_colours,
          show_pct = FALSE,
          show_column_names = FALSE,
          top_annotation = HeatmapAnnotation(
            Age = pat.anno,
            col = list(Age = circlize::colorRamp2(
              c(age_mid - age_half, age_mid, age_mid + age_half),
              c("#0072B2", "white", "#D55E00"))),
            annotation_legend_param = list(Age = list(title = "Age"))),
          row_order = NULL, column_order = colnames(s11.mut.order),
          column_title = NULL,
          heatmap_legend_param = list(title = "Mutation", at = "Mut",
                                      labels = "Mutation"),
          right_annotation = NULL)
dev.off()

## figureS12.tiff — Wilcoxon boxplot: age by mutation status, driver genes

plot_genes <- wilcox_res |> filter(n_mutated >= 3) |> arrange(pval) |> pull(gene)

wx_long <- lapply(plot_genes, function(g) {
  tibble(gene   = g,
         Age    = cli.filt$Age,
         status = if_else(mut.mat.filt[g, ] == 1, "Mutated", "Unmutated"))
}) |> bind_rows() |>
  left_join(wilcox_res |> dplyr::select(gene, pval, padj, n_mutated), by = "gene") |>
  mutate(
    label = sprintf("%s\n(n=%d, p=%.3f)", gene, n_mutated, pval),
    label = fct_reorder(label, pval)
  )

figureS12 <- ggplot(wx_long, aes(status, Age, fill = status)) +
  geom_boxplot(outlier.size = 0.8, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.4) +
  scale_fill_manual(values = c("Mutated" = "#D55E00", "Unmutated" = "#0072B2"),
                    guide = "none") +
  facet_wrap(~ label, nrow = 3, scales = "free_x") +
  labs(x = NULL, y = "Age at diagnosis (years)") +
  my_theme +
  theme(strip.text   = element_text(size = 7),
        axis.text.x  = element_text(size = 7))

ggsave("figureS12.tiff", figureS12,
       width = two_column_width, height = 12, units = "cm", dpi = 300,
       compression = "lzw")
