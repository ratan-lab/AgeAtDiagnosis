
## Date: July 29, 2026 

## RNA-seq age analysis for T-LGL leukemia

# ----------------- Load packages and data --------------------
    new_LIB ="/standard/vol169/cphg_ratan/cphg-RLscratch/share/R/local/4.3.1"
    .libPaths(new_LIB)
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


# set paramters 
  SOFT_POWER      <- 12    # signed network; pickSoftThreshold() confirms R² = 0.945
  MERGE_CUT       <- 0.25
  MIN_MODULE_SIZE <- 30
  N_BOOT          <- 1000
  RERUN_WGCNA     <- FALSE  # set TRUE to rebuild TOM from scratch

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


# ----------------- 1. LGL burden regression + WGCNA   --------------------

    # 1.1  Load counts + clinical, VST normalize 
    rna_path = "/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/RNA-seq"
    count_mat <- read.table(file.path(rna_path, "counts.tsv.gz"))
    rownames(count_mat) <- sub("[.][0-9]+$", "", rownames(count_mat))

    cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/clinical_final.csv")  %>% 
      filter(LGL_Type == "T")
    
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
            TRUE               ~ sample_date >= tx1_date)
      ) |>
        filter(!isTRUE(treated_at_collection)) |>
        column_to_rownames("RegID")

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

    # 1.2 Regress out LGL burden 
      vst_resid <- removeBatchEffect(vst_mat, covariates = meta$cd3cd8)
      input_mat <- t(vst_resid)   # samples × genes

    # 1.3  WGCNA (cached) ─────────────────────────────────────────────────────
    
    wgcna_cache <- file.path(rna_path, "lgl_regress_wgcna.RData") 

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

    # 1.4  Module–trait correlations 
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

    # 1.5  Gene info + identify ribosomal and PCDH-γ modules 
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

    ribo_module <- gene_info |> # annotated as "drakred" by WGCNA 
      filter(grepl("^RPL|^RPS|^FAU$|^UBA52$", symbol)) |>
      dplyr::count(moduleColor, sort = TRUE) |>
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
      dplyr::count(moduleColor, sort = TRUE) |>
      slice_head(n = 1) |>
      pull(moduleColor) # Bisque4 color 
    
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

# --------------------   Figure 3A ---------------------------------
    p_fig3a <- ggplot(eigen_age_df, aes(x = Age, y = eigen)) +
      geom_point(color = "black", alpha = 0.7, size = 1.8) +
      geom_smooth(method = "lm", color = "black", fill = "grey30",
                  alpha = 0.15, linewidth = 0.8) +
      annotate("text", x = Inf, y = Inf,
              label = .ribo_lbl, hjust = 1.1, vjust = 1.5, size = 2.5) +
      labs(x = "Age at diagnosis (years)",
          y = "Ribosomal module eigengene") + 
      my_theme + 
      theme_classic() + 
      theme(axis.title = element_text(size =7))

    ggsave(file.path(fig.path, "Figure3A.png"),p_fig3a, width = 8, height = 6, units = "cm", dpi = 300)


# --------------------   Figure S7. lollipop + compact heatmap  ---------------------------------
    
    cb_cols <- c("Ribosomal" = "#D55E00", "PCDH-γ" = "#009E73", "Other" = "#0072B2")

    mod_df <- tibble(
          module_ME   = rownames(mt_cor),
          module_name = sub("^ME", "", rownames(mt_cor)),
          r_age       = mt_cor[, "Age"],
          p_age       = mt_pval[, "Age"]) |>
      filter(module_name != "grey") |>
      mutate(module_type = factor(case_when(
          module_ME == ribo_ME ~ "Ribosomal",
          module_ME == pcdh_ME ~ "PCDH-γ",
          TRUE                 ~ "Other"
        ), levels = c("Ribosomal", "PCDH-γ", "Other"))) |>
      arrange(r_age) |>
      mutate(module_name = factor(module_name, levels = module_name))

p_lollipop <- ggplot(mod_df, aes(x = r_age, y = module_name)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_segment(aes(x = 0, xend = r_age, yend = module_name), linewidth = 0.5) +
  geom_point(aes(size = -log10(p_age))) +
  geom_text(
    data = filter(mod_df, module_type != "Other"),
    aes(label = sprintf("r=%.2f, p=%.3f", r_age, p_age),
        hjust = ifelse(r_age < 0, -0.5, 1.6)),
    size = 2, show.legend = FALSE
  ) +
  # scale_color_manual(values = cb_cols, name = NULL) +
  scale_size_continuous(name = expression(-log[10](italic(p))), range = c(1.2,3)) +
  coord_cartesian(xlim = c(-0.65, 0.65), clip = "off") +
  labs(x = "Pearson r (eigengene vs Age)", y = NULL) +
  my_theme + 
  theme_bw(base_size = 7) + 
  theme(axis.text.y  = element_text(size = 7),
        legend.title = element_text(
        size = 6,
        face = "plain"
        ),
        legend.text = element_text(size = 7),
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(5, "mm"),

        legend.spacing.y = unit(0.5, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 1))

ggsave(file.path(fig.path, "FigureS7A.png"),p_lollipop, width = 8, height = 10, units = "cm", dpi = 300)

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
# y_colors <- ifelse(paste0("ME", mod_label_order) == ribo_ME, "#D55E00",
#              ifelse(paste0("ME", mod_label_order) == pcdh_ME, "#009E73", "black"))
y_face   <- ifelse(paste0("ME", mod_label_order) %in% c(ribo_ME, pcdh_ME),
                   "bold", "plain")

max_r <- max(heat_long$abs_r, na.rm = TRUE)

heat_long <- heat_long |>
  mutate(
    # 최댓값의 60% 이상이면 흰색 글씨
    text_col = if_else(
      abs_r >= 0.60 * max_r,
      "white",
      "black"
    )
  )


p_heat <- suppressWarnings(
  ggplot(heat_long, aes(x = trait, y = module_name, fill = abs_r)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label, color = text_col), size = 2) +
    # scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
    #                      midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    scale_fill_gradientn(colors = c("white", "grey65", "black"),
      values = c(0, 0.5, 1),
      limits = c(0, max_r),
      breaks = pretty(c(0, max_r), n = 4),
      name = expression("|Pearson r|"),
      oob = scales::squish
    ) + 
    scale_color_identity() +
    labs(x = NULL, y = NULL, 
         subtitle = "Top 6 modules by |r(Age)|; * p < 0.05") +
    my_theme + 
    theme_bw(base_size =7) + 
    theme(axis.text.y     = element_text(face = y_face, size = 7, color = "black"),
          axis.text.x     = element_text(angle = 30, hjust = 1, color ="black"),
          legend.position = "right", 
          legend.text = element_text(size = 7),
        legend.key.height = unit(0.7, "cm"),
        legend.key.width  = unit(0.3, "cm"),
      ) 
    
)

ggsave(file.path(fig.path, "FigureS7B.png"),p_heat, width = 8, height = 10, units = "cm", dpi = 300)


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
m
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

  # figure S8. 
    p_boot <- ggplot(tibble(median_r = boot_medians), aes(median_r)) +
      geom_histogram(bins = 40, fill = "grey85", alpha = 0.5, color = "grey25") +
      geom_vline(xintercept = median(boot_medians), color = "grey25",
                linewidth = 0.5, linetype = "dotted") +
      geom_vline(xintercept = lgll_median, color = "black",
                linewidth = 0.5, linetype = "solid") +
      annotate("text", x = -0.24, y = 76,
              label = "DICE subsample\nmedian", color = "grey25",
              hjust = 1, size = 2) +
      annotate("text", x = lgll_median + 0.005, y = 76,
              label = "LGLL\nmedian", color = "black",
              hjust = 0,  size = 2) +
      labs(x    = "Median age correlation across ribosomal module genes",
          y    = "Count",
          subtitle = sprintf(
              "%d bootstraps (n=%d); DICE median < LGLL: %.1f%%",
              N_BOOT, BOOT_N, 100 * prop_dice_below_lgll)) +
      my_theme + 
      theme_bw(base_size = 7) 

    ggsave(file.path(fig.path, "FigureS8.png"),p_boot, width = 8, height = 6, units = "cm", dpi = 300)


    ## ── 2.5  Trajectory by age bins ──────────────────────────────────────────────────────────
    shared_breaks <- c(17, 29, 39, 49, 59, 69, 79, 89)
    shared_labels <- c("20s", "30s", "40s", "50s", "60s", "70s", "80s")

    # compute the representative values from DICE (for Ribosomal module genes)
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
      mutate(dataset = "LGLL\n(LGL burden adjusted)",
            age_bin = cut(age, breaks = shared_breaks, labels = shared_labels))

    # scale the values within each dataset 
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

    # Figure 3C. 
      p_traj2 <- ggplot(gap_summary,
                        aes(age_bin, mean_z, color = dataset, fill = dataset,
                            group = dataset)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
                    alpha = 0.15, color = NA) +
        geom_line(linewidth = 0.5) +
        geom_point(size = 2, shape = 21) +
        scale_fill_manual(values  = c("DICE CD8+ naive\n(healthy)" = "grey65",
                                      "LGLL\n(LGL burden adjusted)" = "black")) +
        scale_color_manual(values = c("DICE CD8+ naive\n(healthy)" = "grey65",
                                      "LGLL\n(LGL burden adjusted)" = "black")) +
        labs(x = "Age group", y = "Ribosomal module eigengene z-score",
            color = NULL, fill = NULL) +
            my_theme +
        theme_bw(base_size =7) + theme(legend.position = "bottom")

    ggsave(file.path(fig.path, "Figure3C.png"),p_traj2, width = 16, height = 7, units = "cm", dpi = 300)

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

compare_df$new_group = ifelse(compare_df$p_all < 0.05, "p < 0.05", "Not significant") 
p_sens <- ggplot(compare_df, aes(r_primary, r_all, color = new_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey80") +
  geom_point(size = 1, alpha = 0.6) +
  geom_label_repel(
    data = \(d) d |> filter(module_type != "Other" |
                              abs(r_primary - r_all) > 0.15 |
                              abs(r_primary) > 0.25 |
                              abs(r_all)     > 0.25),
    aes(label = module), size = 2, show.legend = FALSE,
    max.overlaps = 20, label.padding = 0.15
  ) +
  scale_color_manual(
    # values = c("Other" = "grey65", "Ribosomal" = "black","PCDH-gamma" = "black"),
    values = c("p < 0.05" = "black", "Not significant" = "grey65"),
    name   = NULL
  ) +
  coord_fixed() +
  labs(x = sprintf("r(ME, Age) : Clinic-only (n=%d)", n_clinic),
       y = sprintf("r(ME, Age) : Clinic + Kit (n=%d)", n_total)) +
  my_theme + 
  theme_bw(base_size =7) +
  theme(legend.position = "bottom") 
  

ggsave(file.path(fig.path, "FigureS9.png"),p_sens, width = 8, height = 7, units = "cm", dpi = 300)


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
