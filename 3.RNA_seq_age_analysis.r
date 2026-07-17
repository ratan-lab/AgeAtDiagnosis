## RNA-seq age analysis for T-LGL leukemia
##
## Produces:
##   figure3.tiff  — Ribosomal eigengene vs Age | DICE violin | trajectory (Fig 3)
##   figureS7.tiff — Module–age lollipop + trait heatmap (Fig S7)
##   figureS8.tiff — Bootstrap histogram (Fig S8)
##   figureS9.tiff — Clinic+Kit sensitivity scatter (Fig S9)

suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(limma)
  library(data.table)
  library(tidyverse)
  library(ggrepel)
  library(annotables)
  library(fgsea)
  library(msigdbr)
  library(patchwork)
  library(readxl)
  library(lubridate)
})

source("funs.R")

## ── Parameters ────────────────────────────────────────────────────────────────
SOFT_POWER      <- 12    # signed network; pickSoftThreshold() confirms R² = 0.945
MERGE_CUT       <- 0.25
MIN_MODULE_SIZE <- 30
N_BOOT          <- 1000
RERUN_WGCNA     <- TRUE  # set TRUE to rebuild TOM from scratch

################################################################################
## PART 1 — LGL BURDEN REGRESSION + WGCNA
################################################################################

## ── 1.1  Load counts + clinical, VST normalize ───────────────────────────────
count_mat <- read.table(rnaseq.path)
rownames(count_mat) <- sub("[.][0-9]+$", "", rownames(count_mat))

cli  <- as.data.frame(fread(clinical.path))
meta <- cli |>
  filter(Clinic_Kit == "Clinic") |>
  mutate(
    Age        = Age_dx,
    STAT3      = STAT3_status,
    cd3cd8     = suppressWarnings(as.numeric(CD3posCD8pos)),
    ALC        = suppressWarnings(as.numeric(ALC)),
    ANC        = suppressWarnings(as.numeric(ANC)),
    RA         = as.integer(RA == "Y"),
    treated    = as.integer(ever_tx == "Y"),
    RegID      = paste0("X", RegID),
    sample_date = lubridate::mdy(Sample_Date),
    tx1_date    = as.Date(suppressWarnings(as.numeric(Tx1_date_raw_clean)),
                          origin = "1899-12-30"),
    treated_at_collection = case_when(
      ever_tx == "N"     ~ FALSE,
      is.na(tx1_date)    ~ NA,
      is.na(sample_date) ~ NA,
      TRUE               ~ sample_date >= tx1_date
    )
  ) |>
  filter(!isTRUE(treated_at_collection)) |>
  column_to_rownames("RegID")

cat(sprintf("Samples after excluding treated-at-collection: %d\n", nrow(meta)))

common <- intersect(rownames(meta), colnames(count_mat))
meta   <- meta[common, ]
counts <- count_mat[, common]

keep   <- !is.na(meta$cd3cd8)
meta   <- meta[keep, ]
counts <- counts[, keep]
cat("Samples:", ncol(counts), "(dropped 1 for missing CD3+CD8+)\n")

BOOT_N <- ncol(counts)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
keep_genes <- rowSums(counts(dds) >= 10) >= ceiling(0.25 * ncol(dds))
dds        <- dds[keep_genes, ]
vst_mat    <- assay(vst(dds, blind = TRUE))
cat("Genes after count filter:", nrow(vst_mat), "\n")

## ── 1.2  Regress out LGL burden ──────────────────────────────────────────────
vst_resid <- removeBatchEffect(vst_mat, covariates = meta$cd3cd8)
input_mat <- t(vst_resid)   # samples × genes

## ── 1.3  WGCNA (cached) ──────────────────────────────────────────────────────
wgcna_cache <- "lgl_regress_wgcna.RData"

if (!RERUN_WGCNA && file.exists(wgcna_cache)) {
  cat("Loading WGCNA from cache:", wgcna_cache, "\n")
  load(wgcna_cache)
} else {
  cat("Building adjacency and TOM (~20 min)...\n")

  sft <- pickSoftThreshold(input_mat, powerVector = seq(6, 20, 2),
                           networkType = "signed", verbose = 3)
  cat(sprintf("Soft power %d: scale-free R² = %.3f (target >= 0.80)\n",
              SOFT_POWER,
              -sign(sft$fitIndices[sft$fitIndices[, 1] == SOFT_POWER, 3]) *
               sft$fitIndices[sft$fitIndices[, 1] == SOFT_POWER, 2]))

  adjacency    <- adjacency(input_mat, power = SOFT_POWER, type = "signed")
  TOM          <- TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM      <- 1 - TOM
  save(dissTOM, file = "lgl_regress_dissTOM.RData")

  gene_tree    <- hclust(as.dist(dissTOM), method = "average")
  dynamic_mods <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = MIN_MODULE_SIZE)
  dynamic_colors <- labels2colors(dynamic_mods)
  merge_res      <- mergeCloseModules(input_mat, dynamic_colors,
                                      cutHeight = MERGE_CUT, verbose = 0)
  module_colors  <- merge_res$colors
  MEs            <- orderMEs(merge_res$newMEs)

  cat("Modules before merge:", length(unique(dynamic_colors)),
      "| After merge:", length(unique(module_colors)), "\n")

  save(MEs, module_colors, gene_tree, file = wgcna_cache)
}

## ── 1.4  Module–trait correlations ───────────────────────────────────────────
traits <- meta |>
  transmute(Age = Age, STAT3 = as.integer(STAT3 == "MT"), RA = RA,
            Treated = treated, ALC = ALC, ANC = ANC) |>
  mutate(across(everything(), as.numeric))
traits  <- traits[, colSums(!is.na(traits)) >= 40, drop = FALSE]
n_obs   <- colSums(!is.na(traits))
mt_cor  <- cor(MEs, traits, use = "pairwise.complete.obs")
mt_pval <- matrix(
  mapply(corPvalueStudent, as.vector(mt_cor), rep(n_obs, each = nrow(mt_cor))),
  nrow = nrow(mt_cor), dimnames = dimnames(mt_cor)
)
mt_padj <- matrix(p.adjust(mt_pval, method = "BH"),
                  nrow = nrow(mt_pval), dimnames = dimnames(mt_pval))

## ── 1.5  Gene info + identify ribosomal and PCDH-γ modules ──────────────────
modNames <- sub("^ME", "", names(MEs))
mm_df    <- as.data.frame(cor(input_mat, MEs, use = "p"))
names(mm_df) <- paste0("MM.", modNames)

meta_aligned <- meta[rownames(input_mat), ]
gene_info <- tibble(
  ensGene     = colnames(input_mat),
  moduleColor = module_colors,
  GS_age      = cor(input_mat, meta_aligned$Age)[, 1]
) |>
  left_join(grch38 |> distinct(ensgene, .keep_all = TRUE) |> select(ensgene, symbol),
            by = c("ensGene" = "ensgene")) |>
  filter(!is.na(symbol), symbol != "") |>
  distinct(ensGene, .keep_all = TRUE)

gene_info <- bind_cols(gene_info,
                       mm_df[match(gene_info$ensGene, colnames(input_mat)), ])

ribo_module <- gene_info |>
  filter(grepl("^RPL|^RPS|^FAU$|^UBA52$", symbol)) |>
  count(moduleColor, sort = TRUE) |>
  slice_head(n = 1) |>
  pull(moduleColor)
stopifnot("no RP genes found in any module" = length(ribo_module) == 1L)
ribo_ME <- paste0("ME", ribo_module)

cat(sprintf("\nRibosomal module: %s (%d RP genes)  r(Age)=%.3f  p=%.3f\n",
            ribo_module,
            sum(grepl("^RPL|^RPS|^FAU$|^UBA52$",
                      gene_info$symbol[gene_info$moduleColor == ribo_module])),
            mt_cor[ribo_ME, "Age"], mt_pval[ribo_ME, "Age"]))

orange_genes <- gene_info |>
  filter(moduleColor == ribo_module) |>
  select(ensGene, symbol, GS_age)

pcdh_module <- gene_info |>
  filter(grepl("^PCDHG", symbol)) |>
  count(moduleColor, sort = TRUE) |>
  slice_head(n = 1) |>
  pull(moduleColor)
stopifnot("no PCDHG genes found in any module" = length(pcdh_module) == 1L)
pcdh_ME <- paste0("ME", pcdh_module)

cat(sprintf("PCDH-γ module: %s (%d PCDHG genes)  r(Age)=%.3f  p=%.3f\n",
            pcdh_module,
            sum(grepl("^PCDHG", gene_info$symbol[gene_info$moduleColor == pcdh_module])),
            mt_cor[pcdh_ME, "Age"], mt_pval[pcdh_ME, "Age"]))

## ── 1.5b  Figure 3A: ribosomal eigengene vs Age scatter ─────────────────────
stopifnot("MEs rows not all in meta" = all(rownames(MEs) %in% rownames(meta)))
eigen_age_df <- tibble(
  Age   = meta[rownames(MEs), "Age"],
  eigen = MEs[, ribo_ME]
)

.ribo_pv  <- mt_pval[ribo_ME, "Age"]
.ribo_lbl <- sprintf("r = %.3f\n%s", mt_cor[ribo_ME, "Age"],
                     if (.ribo_pv < 0.001) "p < 0.001"
                     else sprintf("p = %.3f", .ribo_pv))

p_fig3a <- ggplot(eigen_age_df, aes(x = Age, y = eigen)) +
  geom_point(color = "#D55E00", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#D55E00", fill = "#D55E00",
              alpha = 0.15, linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = .ribo_lbl, hjust = 1.1, vjust = 1.5, size = 3) +
  labs(x = "Age at diagnosis (years)",
       y = "Ribosomal module eigengene") +
  theme_bw(base_size = 10)

## ── 1.6  Figure S7: lollipop + compact heatmap ───────────────────────────────
cb_cols <- c("Ribosomal" = "#D55E00", "PCDH-γ" = "#009E73", "Other" = "#0072B2")

mod_df <- tibble(
  module_ME   = rownames(mt_cor),
  module_name = sub("^ME", "", rownames(mt_cor)),
  r_age       = mt_cor[, "Age"],
  p_age       = mt_pval[, "Age"]
) |>
  filter(module_name != "grey") |>
  mutate(module_type = factor(case_when(
    module_ME == ribo_ME ~ "Ribosomal",
    module_ME == pcdh_ME ~ "PCDH-γ",
    TRUE                 ~ "Other"
  ), levels = c("Ribosomal", "PCDH-γ", "Other"))) |>
  arrange(r_age) |>
  mutate(module_name = factor(module_name, levels = module_name))

p_lollipop <- ggplot(mod_df, aes(x = r_age, y = module_name, color = module_type)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_segment(aes(x = 0, xend = r_age, yend = module_name), linewidth = 0.5) +
  geom_point(aes(size = -log10(p_age))) +
  geom_text(
    data = filter(mod_df, module_type != "Other"),
    aes(label = sprintf("r=%.2f, p=%.3f", r_age, p_age),
        hjust = ifelse(r_age < 0, 1.12, -0.12)),
    size = 2.8, show.legend = FALSE
  ) +
  scale_color_manual(values = cb_cols, name = NULL) +
  scale_size_continuous(name = expression(-log[10](italic(p))), range = c(2, 5)) +
  coord_cartesian(xlim = c(-0.65, 0.65), clip = "off") +
  labs(x = "Pearson r (eigengene vs Age)", y = NULL, tag = "A") +
  theme_bw(base_size = 10) +
  theme(axis.text.y  = element_text(size = 7),
        plot.tag     = element_text(face = "bold"),
        plot.margin  = margin(5, 60, 5, 5))

top_mods <- mod_df |>
  arrange(desc(abs(r_age))) |>
  slice_head(n = 6) |>
  arrange(r_age) |>
  pull(module_ME) |>
  as.character()

heat_long <- as_tibble(mt_cor[top_mods, , drop = FALSE], rownames = "module_ME") |>
  pivot_longer(-module_ME, names_to = "trait", values_to = "r") |>
  left_join(
    as_tibble(mt_pval[top_mods, , drop = FALSE], rownames = "module_ME") |>
      pivot_longer(-module_ME, names_to = "trait", values_to = "p"),
    by = c("module_ME", "trait")
  ) |>
  mutate(
    label       = ifelse(p < 0.05, sprintf("%.2f*", r), sprintf("%.2f", r)),
    module_name = sub("^ME", "", module_ME),
    module_type = factor(case_when(
      module_ME == ribo_ME ~ "Ribosomal",
      module_ME == pcdh_ME ~ "PCDH-γ",
      TRUE                 ~ "Other"
    ), levels = c("Ribosomal", "PCDH-γ", "Other")),
    module_name = factor(module_name, levels = sub("^ME", "", top_mods)),
    trait       = factor(trait, levels = colnames(mt_cor))
  )

mod_label_order <- levels(heat_long$module_name)
y_colors <- ifelse(paste0("ME", mod_label_order) == ribo_ME, "#D55E00",
             ifelse(paste0("ME", mod_label_order) == pcdh_ME, "#009E73", "black"))
y_face   <- ifelse(paste0("ME", mod_label_order) %in% c(ribo_ME, pcdh_ME),
                   "bold", "plain")

p_heat <- suppressWarnings(
  ggplot(heat_long, aes(x = trait, y = module_name, fill = r)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.8) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    labs(x = NULL, y = NULL, tag = "B",
         subtitle = "Top 6 modules by |r(Age)|; * p < 0.05") +
    theme_bw(base_size = 10) +
    theme(axis.text.y     = element_text(color = y_colors, face = y_face, size = 8),
          axis.text.x     = element_text(angle = 30, hjust = 1),
          plot.tag        = element_text(face = "bold"),
          legend.position = "right")
)

## ── ORA: ribosomal module pathway enrichment ──────────────────────────────────
h_sets    <- msigdbr(species = "human", collection = "H")
k_sets    <- msigdbr(species = "human", collection = "C2",
                     subcollection = "CP:KEGG_LEGACY")
msig_list <- c(split(h_sets$gene_symbol, h_sets$gs_name),
               split(k_sets$gene_symbol, k_sets$gs_name))

orange_set <- unique(orange_genes$symbol)
universe   <- unique(gene_info$symbol)

ora_res <- fora(
  pathways = msig_list,
  genes    = orange_set,
  universe = universe,
  minSize  = 15, maxSize  = 500
) |>
  as_tibble() |>
  filter(padj < 0.05) |>
  mutate(
    pathway      = str_remove(pathway, "^HALLMARK_|^KEGG_"),
    enrich_ratio = (overlap / length(orange_set)) / (size / length(universe))
  ) |>
  arrange(desc(enrich_ratio))

cat(sprintf("ORA: %d pathways significant (FDR<0.05) in ribosomal module\n", nrow(ora_res)))
print(ora_res |> select(pathway, overlap, size, padj, enrich_ratio) |> head(10))

################################################################################
## PART 2 — DICE CD8+ NAIVE COMPARISON + BOOTSTRAP
################################################################################

## ── 2.1  DICE metadata and TPM ───────────────────────────────────────────────
meta_dice <- read_xlsx(dice.metadata.path, skip = 2) |>
  filter(!is.na(`Donor ID`), !grepl("_2$", `Donor ID`)) |>
  mutate(donor_id = as.integer(`Donor ID`), age = as.numeric(`Age (y)`)) |>
  filter(!is.na(age), donor_id <= 89) |>
  mutate(col_idx = as.character(donor_id - 1)) |>
  select(col_idx, donor_id, age)

cat("\nDICE donors:", nrow(meta_dice),
    "| Age range:", min(meta_dice$age), "-", max(meta_dice$age), "\n")

dice_header     <- fread(dice.cd8.path, nrows = 0)
dice_avail_cols <- setdiff(names(dice_header), "Feature_name")
dice_missing    <- setdiff(meta_dice$col_idx, dice_avail_cols)
if (length(dice_missing) > 0)
  stop(sprintf("DICE col_idx mismatch for donors: %s",
               paste(dice_missing, collapse = ", ")))

tpm_raw <- fread(dice.cd8.path, select = c("Feature_name", meta_dice$col_idx)) |>
  mutate(ensGene = sub("[.][0-9]+$", "", Feature_name)) |>
  filter(ensGene %in% orange_genes$ensGene) |>
  select(-Feature_name)

stopifnot(all(meta_dice$col_idx %in% names(tpm_raw)))

## ── 2.2  Per-gene age correlations + Wilcoxon ────────────────────────────────
count_dice <- tpm_raw |> column_to_rownames("ensGene") |>
  select(all_of(meta_dice$col_idx)) |> as.matrix()
log_tpm  <- log2(count_dice + 1)
age_dice <- meta_dice$age[match(colnames(log_tpm), meta_dice$col_idx)]

suppressWarnings(
  r_dice_raw <- cor(t(log_tpm), age_dice, use = "pairwise.complete.obs")[, 1]
)

dice_cors <- tibble(ensGene = names(r_dice_raw), r_dice = r_dice_raw) |>
  filter(!is.na(r_dice)) |>
  inner_join(orange_genes, by = "ensGene") |>
  mutate(is_RP = grepl("^RPL|^RPS|^FAU$|^UBA52$", symbol))

wt_dice <- wilcox.test(dice_cors$r_dice, dice_cors$GS_age, paired = TRUE)
cat(sprintf("\nDICE naive median r = %.3f  |  LGLL GS_age median = %.3f\n",
            median(dice_cors$r_dice), median(dice_cors$GS_age)))
cat(sprintf("Paired Wilcoxon p = %.2e\n", wt_dice$p.value))

## ── 2.3  Bootstrap ───────────────────────────────────────────────────────────
set.seed(42)
boot_idx     <- match(dice_cors$ensGene, rownames(log_tpm))
stopifnot(sum(is.na(boot_idx)) == 0)
boot_mat     <- log_tpm[boot_idx, ]
boot_age     <- meta_dice$age[match(colnames(boot_mat), meta_dice$col_idx)]
lgll_median  <- median(dice_cors$GS_age, na.rm = TRUE)

boot_medians <- map_dbl(seq_len(N_BOOT), function(i) {
  idx   <- sample(ncol(boot_mat), BOOT_N, replace = FALSE)
  r_sub <- cor(t(boot_mat[, idx]), boot_age[idx],
               use = "pairwise.complete.obs")[, 1]
  median(r_sub, na.rm = TRUE)
})

prop_dice_below_lgll <- mean(boot_medians < lgll_median)
cat(sprintf(
  "\nBootstrap (%d iterations, n=%d): DICE subsample median = %.3f [%.3f, %.3f]\n",
  N_BOOT, BOOT_N, median(boot_medians),
  quantile(boot_medians, 0.025), quantile(boot_medians, 0.975)))
cat(sprintf("%.1f%% of bootstraps: DICE subsample median < LGLL median\n",
            100 * prop_dice_below_lgll))

## ── 2.4  Plots ───────────────────────────────────────────────────────────────
p_dice_violin <- dice_cors |>
  select(symbol, is_RP, r_dice, GS_age) |>
  pivot_longer(c(r_dice, GS_age), names_to = "dataset", values_to = "r_age") |>
  mutate(dataset = case_when(
    dataset == "r_dice" ~ "DICE CD8+ naive\n(healthy)",
    dataset == "GS_age" ~ "LGLL\n(LGL burden regressed)"
  )) |>
  ggplot(aes(dataset, r_age, fill = dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = c("DICE CD8+ naive\n(healthy)" = "#0072B2",
                                "LGLL\n(LGL burden regressed)" = "#D55E00")) +
  labs(x = NULL, y = "Pearson r with age") +
  theme_bw() + theme(legend.position = "none")

p_boot <- ggplot(tibble(median_r = boot_medians), aes(median_r)) +
  geom_histogram(bins = 40, fill = "#0072B2", alpha = 0.7, color = "white") +
  geom_vline(xintercept = median(boot_medians), color = "#0072B2",
             linewidth = 1, linetype = "solid") +
  geom_vline(xintercept = lgll_median, color = "#D55E00",
             linewidth = 1, linetype = "dashed") +
  annotate("text", x = median(boot_medians) - 0.003, y = Inf,
           label = "DICE\nsubsample\nmedian", color = "#0072B2",
           hjust = 1, vjust = 1.5, size = 3) +
  annotate("text", x = lgll_median + 0.003, y = Inf,
           label = "LGLL\nmedian", color = "#D55E00",
           hjust = 0, vjust = 1.5, size = 3) +
  labs(x    = "Median age correlation across ribosomal module genes",
       y    = "Count",
       subtitle = sprintf(
         "%d iterations, n=%d subsampled  |  %.1f%% bootstraps: DICE median < LGLL median",
         N_BOOT, BOOT_N, 100 * prop_dice_below_lgll)) +
  theme_bw()

## ── 2.5  Trajectory ──────────────────────────────────────────────────────────
shared_breaks <- c(17, 29, 39, 49, 59, 69, 79, 89)
shared_labels <- c("20s", "30s", "40s", "50s", "60s", "70s", "80s")

pca_dice   <- prcomp(t(log_tpm[rowSums(log_tpm) > 0, ]), center = TRUE, scale. = TRUE)
dice_eigen <- pca_dice$x[, 1]
if (cor(dice_eigen, colMeans(log_tpm)) < 0) dice_eigen <- -dice_eigen

lgll_ages_traj <- meta[rownames(MEs), "Age"]

dice_traj <- tibble(col_idx = names(dice_eigen), eigen = dice_eigen) |>
  left_join(meta_dice, by = "col_idx") |>
  mutate(dataset = "DICE CD8+ naive\n(healthy)",
         age_bin = cut(age, breaks = shared_breaks, labels = shared_labels))
lgll_traj <- tibble(col_idx = rownames(MEs), eigen = MEs[, ribo_ME],
                    age = lgll_ages_traj) |>
  mutate(dataset = "LGLL\n(LGL burden regressed)",
         age_bin = cut(age, breaks = shared_breaks, labels = shared_labels))

traj_df <- bind_rows(
  dice_traj |> mutate(eigen_z = scale(eigen)[, 1]) |> select(dataset, age_bin, eigen_z),
  lgll_traj |> mutate(eigen_z = scale(eigen)[, 1]) |> select(dataset, age_bin, eigen_z)
)

gap_summary <- traj_df |>
  filter(!is.na(age_bin)) |>
  group_by(dataset, age_bin) |>
  summarise(mean_z = mean(eigen_z, na.rm = TRUE),
            se_z   = sd(eigen_z,   na.rm = TRUE) / sqrt(n()),
            n      = n(), .groups = "drop")

p_traj2 <- ggplot(gap_summary,
                   aes(age_bin, mean_z, color = dataset, fill = dataset,
                       group = dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values  = c("DICE CD8+ naive\n(healthy)" = "#0072B2",
                                 "LGLL\n(LGL burden regressed)" = "#D55E00")) +
  scale_color_manual(values = c("DICE CD8+ naive\n(healthy)" = "#0072B2",
                                 "LGLL\n(LGL burden regressed)" = "#D55E00")) +
  labs(x = "Age group", y = "Ribosomal module eigengene (z-score, mean ± SE)",
       color = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "bottom")

################################################################################
## PART 3 — SENSITIVITY: Clinic + Kit samples
################################################################################

cat("\n=== PART 3: Sensitivity — Clinic + Kit samples ===\n")

## ── 3.1  All-sample metadata ─────────────────────────────────────────────────
meta_all <- cli |>
  filter(Clinic_Kit %in% c("Clinic", "Kit")) |>
  mutate(
    Age        = Age_dx,
    STAT3      = STAT3_status,
    cd3cd8     = suppressWarnings(as.numeric(CD3posCD8pos)),
    ALC        = suppressWarnings(as.numeric(ALC)),
    ANC        = suppressWarnings(as.numeric(ANC)),
    RA         = as.integer(RA == "Y"),
    treated    = as.integer(ever_tx == "Y"),
    clinic_kit = as.integer(Clinic_Kit == "Kit"),
    RegID      = paste0("X", RegID)
  ) |>
  column_to_rownames("RegID") |>
  filter(!is.na(cd3cd8))

common_all <- intersect(rownames(meta_all), colnames(count_mat))
meta_all   <- meta_all[common_all, ]
counts_all <- count_mat[, common_all]

cat(sprintf("All-sample cohort: n=%d  (Clinic=%d, Kit=%d)\n",
            nrow(meta_all),
            sum(meta_all$clinic_kit == 0),
            sum(meta_all$clinic_kit == 1)))

## ── 3.2  VST + regress out cd3cd8 and clinic_kit ────────────────────────────
dds_all <- DESeqDataSetFromMatrix(countData = counts_all,
                                   colData   = meta_all,
                                   design    = ~ 1)
dds_all     <- dds_all[rowSums(counts(dds_all) >= 10) >= ceiling(0.25 * ncol(dds_all)), ]
vst_all     <- assay(vst(dds_all, blind = TRUE))
vst_all_res <- removeBatchEffect(
  vst_all,
  covariates = cbind(meta_all$cd3cd8, meta_all$clinic_kit)
)
input_all <- t(vst_all_res)

## ── 3.3  Project module assignments onto all-sample matrix ───────────────────
if (is.null(names(module_colors))) names(module_colors) <- colnames(input_mat)
genes_shared <- intersect(colnames(input_all), names(module_colors))
cat(sprintf("Shared genes: %d / %d\n", length(genes_shared), length(module_colors)))

MEs_all <- moduleEigengenes(input_all[, genes_shared],
                             colors    = module_colors[genes_shared],
                             softPower = SOFT_POWER)$eigengenes
MEs_all <- orderMEs(MEs_all)
stopifnot(all(c(ribo_ME, pcdh_ME) %in% names(MEs_all)))

## ── 3.4  Module–Age correlations in expanded cohort ─────────────────────────
age_all <- meta_all[rownames(MEs_all), "Age"]
n_age   <- sum(!is.na(age_all))
r_all   <- cor(MEs_all, age_all, use = "pairwise.complete.obs")[, 1]
p_all   <- sapply(r_all, function(r) corPvalueStudent(r, n_age))

cat(sprintf("Ribosomal  — Primary: r=%.3f p=%.3f | Sensitivity: r=%.3f p=%.3f\n",
            mt_cor[ribo_ME, "Age"], mt_pval[ribo_ME, "Age"],
            r_all[ribo_ME], p_all[ribo_ME]))
cat(sprintf("PCDH-gamma — Primary: r=%.3f p=%.3f | Sensitivity: r=%.3f p=%.3f\n",
            mt_cor[pcdh_ME, "Age"], mt_pval[pcdh_ME, "Age"],
            r_all[pcdh_ME], p_all[pcdh_ME]))

## ── 3.5  Sensitivity scatter ─────────────────────────────────────────────────
shared_mods <- intersect(sub("^ME", "", rownames(mt_cor)),
                          sub("^ME", "", names(MEs_all)))
compare_df <- tibble(
  module      = shared_mods,
  r_primary   = mt_cor [paste0("ME", shared_mods), "Age"],
  p_primary   = mt_pval[paste0("ME", shared_mods), "Age"],
  r_all       = r_all  [paste0("ME", shared_mods)],
  p_all       = p_all  [paste0("ME", shared_mods)],
  module_type = factor(case_when(
    shared_mods == ribo_module ~ "Ribosomal",
    shared_mods == pcdh_module ~ "PCDH-gamma",
    TRUE                       ~ "Other"
  ), levels = c("Ribosomal", "PCDH-gamma", "Other"))
)

n_clinic <- sum(meta_all$clinic_kit == 0)
n_total  <- nrow(meta_all)

p_sens <- ggplot(compare_df, aes(r_primary, r_all, color = module_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_point(size = 3, alpha = 0.85) +
  geom_label_repel(
    data = \(d) d |> filter(module_type != "Other" |
                              abs(r_primary - r_all) > 0.15 |
                              abs(r_primary) > 0.25 |
                              abs(r_all)     > 0.25),
    aes(label = module), size = 2.8, show.legend = FALSE,
    max.overlaps = 20, label.padding = 0.15
  ) +
  scale_color_manual(
    values = c("Other" = "#0072B2", "Ribosomal" = "#D55E00",
               "PCDH-gamma" = "#009E73"),
    name   = NULL
  ) +
  coord_fixed() +
  labs(x = sprintf("r(ME, Age) — Clinic-only (n=%d)", n_clinic),
       y = sprintf("r(ME, Age) — Clinic + Kit (n=%d)", n_total)) +
  theme_bw() +
  theme(legend.position = "bottom")

################################################################################
## MANUSCRIPT FIGURES
################################################################################

p3a <- p_fig3a       + my_theme
p3b <- p_dice_violin + my_theme + theme(legend.position = "none")
p3c <- p_traj2       + my_theme + theme(legend.position = "bottom")

ggsave("figure3.tiff",
       (p3a | p3b) / p3c + plot_annotation(tag_levels = "A"),
       width = 17, height = 14, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figure3.tiff\n")

p_s7_ms <- (p_lollipop + my_theme) | suppressWarnings(p_heat + my_theme)
ggsave("figureS7.tiff", p_s7_ms,
       width = 18, height = 10, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS7.tiff\n")

ggsave("figureS8.tiff", p_boot + my_theme,
       width = onepointfive_column_width, height = 7, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS8.tiff\n")

ggsave("figureS9.tiff",
       p_sens + my_theme + theme(legend.position = "bottom"),
       width = 9, height = 9, units = "cm", dpi = 300, compression = "lzw")
cat("Saved figureS9.tiff\n")
