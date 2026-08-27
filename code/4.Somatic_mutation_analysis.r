## Last updated: July 31,

## This script is for somatic mutation analysis and creating the plots



# ----------------- Load packages and data --------------------
    new_LIB ="/standard/vol169/cphg_ratan/cphg-RLscratch/share/R/local/4.3.1"
    .libPaths(new_LIB)
    library(ComplexHeatmap)
    library(tidyverse)
    library(data.table)
    library(logistf)

    source("funs.R")

    ## ── Somatic mutation data ───────────────────────────────────────────────────
    dat.path = "/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data"
    somatic.path <- file.path(dat.path, "mutation", "somatic_mutations.tsv.gz") 
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


    ## ── Clinical data ───────────────────────────────────────────────────
    # cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/Age_heejin_Nate_v4_07102026.csv") 
    cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/clinical.csv") 
    
    cli = cli %>% 
        rename(Age = Age_dx) %>%
        filter(!is.na(patient_id) & patient_id != "") %>%
        mutate(RegID.2 = paste0(patient_id, "_tumor"))

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

    cor_label <- sprintf("Spearman rho = %.2f, p = %.3g\nAge adj",
                     cor_tmb$estimate, cor_tmb$p.value, lm_age_p)

    lm_adj    <- lm(TMB ~ Age + STAT3_status, data = filter(tmb_df, !is.na(STAT3_status)))
    lm_coefs  <- coef(summary(lm_adj))
    lm_age_p  <- lm_coefs["Age", "Pr(>|t|)"]
  
      # Figure 4A. TMB ~ Age 
      figure4a <- ggplot(filter(tmb_df, !is.na(STAT3_status)),aes(Age, TMB))+
          geom_point(alpha = 0.7, color = "black") +
          geom_smooth(method = "lm", se = TRUE, alpha = 0.15, color = "black") +
          # scale_color_manual(values = c(WT = "black", MT = "#D55E00"), name = "STAT3") +
          # scale_fill_manual(values  = c(WT = "black", MT = "#D55E00"), guide = "none") +
          annotate("text", x = 18, y = 180, label = cor_label,
                  hjust = 0, vjust = 1.5, size = 2) +
          labs(x = "Age at diagnosis (years)",
              y = "Tumor mutation burden (TMB)") +
          my_theme + 
          theme_classic(base_size = 7) 


      fig.path ="//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/figures"
      ggsave(file.path(fig.path, "Figure4A.png"),figure4a, width = 8, height = 6, units = "cm", dpi = 300)


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
      sig.genes <- df.genes |> filter(padj < 0.05) |> pull(gene) # nothing 


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

    ## figure4.tiff — TMB vs age, colored by STAT3 status

    cor_label <- sprintf("Spearman rho = %.2f, p = %.3g Age adj. \nfor STAT3: p = %.3g",
                        cor_tmb$estimate, cor_tmb$p.value, lm_age_p)

    figure4b<- ggplot(filter(tmb_df, !is.na(STAT3_status)),
                      aes(Age, TMB, color = STAT3_status, fill = STAT3_status)) +
        geom_point(alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
        scale_color_manual(values = c(WT = "grey55", MT = "black"), name = "STAT3") +
        scale_fill_manual(values  = c(WT = "grey55", MT = "black"), guide = "none") +
        annotate("text", x = 15, y= 185, label = cor_label,
              hjust = 0, vjust = 1.5, size = 2) +
      labs(x = "Age at diagnosis (years)",
          y = "Tumor mutation burden (TMB)") +
      my_theme + 
      theme_classic(base_size = 7)

    fig.path ="//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/figures"
    ggsave(file.path(fig.path, "Figure4B.png"),figure4b, width = 8, height = 6, units = "cm", dpi = 300)


## figureS10.tiff — Driver gene OncoPrint sorted by age (14 x 6 cm)
##   restricted to LGLL driver genes with >= 5 mutated patients

      s10.genes <- intersect(TLGL_DRIVER_GENES, names(which(rowSums(mut.mat.filt) >= 5)))
      s10.mut   <- mut.mat.filt[s10.genes, ]

      s10.mut <- matrix(
        ifelse(s10.mut == 1, "Mut", " "),
        nrow = nrow(s10.mut), ncol = ncol(s10.mut),
        dimnames = dimnames(s10.mut)
      )

      mut_colours <- c("Mut" = "black")
      alter_fun <- list(
        background = function(x, y, w, h) {
          grid.rect(x, y, 
            # w - unit(0.5, "mm"), h - unit(0.5, "mm"),
            width= w*0.4, height = h *0.95, 
                    gp = gpar(fill = "grey85", col = NA))
        },
        Mut = function(x, y, w, h) {
          grid.rect(x, y, 
            # w - unit(0.5, "mm"), h - unit(0.5, "mm"),
            width= w*0.4, height = h *0.95,
                    gp = gpar(fill = mut_colours["Mut"], col = NA))
        }
      )

      pat.anno.dat <- data.frame(Age = cli.filt$Age, ID = cli.filt$RegID.2)
      pat.anno.dat <- pat.anno.dat[order(pat.anno.dat$Age), ]
      pat.anno <- pat.anno.dat$Age
      names(pat.anno) <- pat.anno.dat$ID

      col_order <- match(names(pat.anno), colnames(s10.mut))
      stopifnot(!anyNA(col_order))
      s10.mut.order <- s10.mut[, col_order]

      age_min = min(pat.anno) 
      age_mid  <- median(pat.anno)
      age_max = max(pat.anno)
      # age_half <- max(abs(pat.anno - age_mid))

      png(file.path(fig.path, "FigureS10.png"), width = 8, height =4, units ="cm", res =300)
      oncoPrint(s10.mut.order,
                alter_fun = alter_fun, col = mut_colours,
                show_pct = FALSE,
                show_column_names = FALSE,
                row_names_gp = grid::gpar(fontsize = 6), 
                
                top_annotation = HeatmapAnnotation(
                  Age = pat.anno,
                  col = list(Age = circlize::colorRamp2(
                    # c(age_mid - age_half, age_mid, age_mid + age_half),
                    c(age_min, age_mid, age_max), 
                    c("white", "grey70", "black"))),
                    annotation_legend_param = list(
                          Age = list(
                            title = "Age",
                            at = c(age_min, round(age_mid), age_max),
                            labels = c(as.character(age_min), as.character(round(age_mid)), as.character(age_max)),
                            title_gp  = grid::gpar(fontsize = 6),
                            labels_gp = grid::gpar(fontsize = 5),
                            grid_width  = grid::unit(4, "mm"),
                            grid_height = grid::unit(2, "mm"),
                            title_gap = grid::unit(1, "mm")
                          )
                        ),
                        annotation_name_gp = grid::gpar(fontsize = 6)
                      ),

                row_order = NULL, column_order = colnames(s10.mut.order),
                column_title = NULL,
                  heatmap_legend_param = list(
                          title = "Mutation",
                          at = "Mut",
                          labels = "Mutation",
                          title_gp  = grid::gpar(fontsize = 6),
                          labels_gp = grid::gpar(fontsize = 5),
                          grid_width  = grid::unit(2, "mm"),
                          grid_height = grid::unit(2, "mm"),
                          title_gap = grid::unit(1, "mm")
                        ),
                right_annotation = NULL)
      dev.off()

## figureS11.tiff — Wilcoxon boxplot: age by mutation status, driver genes

plot_genes <- wilcox_res |> 
        # filter(n_mutated >= 3) |> 
        arrange(pval) |> pull(gene)

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

figureS11 <- ggplot(wx_long, aes(status, Age, fill = status)) +
  geom_boxplot(outlier.size = 0.8, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.4) +
  scale_fill_manual(values = c("Mutated" = "black", "Unmutated" = "grey65"),
                    guide = "none") +
  facet_wrap(~ label, nrow = 3, scales = "fixed",axes = "margins", axis.labels = "margins") +
  labs(x = NULL, y = "Age at diagnosis (years)", title = "15 LGLL driver genes") +
  my_theme + 
  theme_bw(base_size = 7) + 
  theme(strip.text   = element_text(size = 7),
        axis.text.x  = element_text(size = 7, angle = 30, hjust = 1), 
        strip.background = element_rect(fill = "white", color ="black", linewidth =0.3), 
        plot.title = element_text(size=8, hjust = 0.5)) 
ggsave(file.path(fig.path, "FigureS11.png"),figureS11, width = 16, height = 12, units = "cm", dpi = 300)
