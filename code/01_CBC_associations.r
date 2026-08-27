

## Date: Ju;y 24, 2026 
## writer: Jisu Shin 


## LGL Leukemia Clinical data: Heejin + Nate Emma data first 


# ----------------- Load packages and data --------------------
library(data.table)
library(tidyverse)
library(lubridate)
library(ggpubr) 
library(haven)
library(cowplot) 
library(cocor) 
library(readxl)
library(patchwork)


# cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/Age_heejin_Nate_v4_07102026.csv") 
cli = read_csv("/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/clinical_final.csv") 
cli = cli %>% filter(LGL_Type == "T")

# ----------------- 0. functions  --------------------

RINT <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

rm.outliers <- function(df) {
  q   <- quantile(df[[1]], probs = c(.25, .75))
  iqr <- IQR(df[[1]])
  df[df[[1]] > q[1] - 1.5 * iqr & df[[1]] < q[2] + 1.5 * iqr, ]
}

point.keepers <- function(df) {
  q   <- quantile(df[[1]], probs = c(.25, .75))
  iqr <- IQR(df[[1]])
  df[[1]] > q[1] - 1.5 * iqr & df[[1]] < q[2] + 1.5 * iqr
}


my_theme <- theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    plot.title      = element_text(size = 8),
    axis.title      = element_text(size = 8),
    axis.text       = element_text(size = 7),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    strip.text      = element_text(size = 8),
    text            = element_text(family = "Arial")
  )

single_column_width <- 8
onepointfive_column_width <- 12
two_column_width <- 16figure4a <- ggplot(filter(tmb_df, !is.na(STAT3_status)),aes(Age, TMB))+
          geom_point(alpha = 0.7, color = "black") +
          geom_smooth(method = "lm", se = TRUE, alpha = 0.15, color = "black") +
          # scale_color_manual(values = c(WT = "black", MT = "#D55E00"), name = "STAT3") +
          # scale_fill_manual(values  = c(WT = "black", MT = "#D55E00"), guide = "none") +
          annotate("text", x = Inf, y = Inf, label = cor_label,
                  hjust = 1.1, vjust = 1.5, size = 2) +
          labs(x = "Age at diagnosis (years)",
              y = "Tumor mutation burden") +
          my_theme

fig.path ="//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/figures"

# ----------------- 0. prep data   --------------------
    
    # T-LGLL 
        # data = read_csv("//standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/manuscript/data/clinical.csv")
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

        cat("T-LGL patients (CBC and treatment QC passed):", nrow(cli.new), "\n")

        patient = cli.new %>% dplyr::select(where(is.numeric), STAT3) %>% 
                    dplyr::select(-RegID)


    # NHANES data 
        nhanes.path = "/standard/vol169/cphg_ratan/cphg-RLscratch/dzs5cf/LGL_age_onset/data/clinical/NHANES"
        nhanes.cbc.path <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/CBC_L.xpt"
        nhanes.demo.path <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.xpt"
    
        nhanes.blood <- read_xpt(nhanes.cbc.path)
        nhanes.demo  <- read_xpt(nhanes.demo.path)
        nhanes.names <- read_csv(file.path(nhanes.path, "CDC_tempus_cbc_corresponding_columns.csv"))

        normal.blood <- nhanes.blood |> dplyr::select(SEQN, all_of(nhanes.names$CDC))
        colnames(normal.blood) <- c("SEQN", nhanes.names$TEMPUS_CBC)

        normal.demo <- nhanes.demo |>
            dplyr::select(SEQN, Age = RIDAGEYR, Sex = RIAGENDR) |>
                mutate(Sex = ifelse(Sex == 1, "M", "F"))

        normal <- left_join(normal.blood, normal.demo, by = "SEQN") |>
            mutate(Group = "Normal")

        cat("NHANES reference N =", nrow(normal), "\n")

        ## CBC variables to analyse.
        vars <- sort(colnames(normal)[c(2, seq(5, (ncol(normal) - 3)))])



# ----------------- Figure S1. associations in T-LGLL    --------------------
## Figure S1 

patient.plots <- list()
patients_res  <- tibble()

for (var in vars) {
  total <- patient |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
  if (nrow(total) < 15) next
  total2 <- total[point.keepers(total), ]
  if (nrow(total2) < 15) next
  total2[[var]] <- RINT(total2[[var]])

  model   <- lm(as.formula(paste0(var, " ~ Age")), data = total2)
  total.p <- coef(summary(model))[2, 4]
  total.b <- coef(summary(model))[2, 1]
  patients_res <- bind_rows(patients_res, tibble(name = var, pval = total.p, beta = total.b))

  p <- ggscatter(total2, x = "Age", y = var, add = "reg.line", conf.int = TRUE,
                  cor.coeff.args = list(method = "pearson"), alpha = 0.3) +
    coord_cartesian(xlim = c(20, 80)) +
    stat_cor(label.x = 30) +
    ggtitle(var) + 
    labs(x = "Age at diagnosis", y = "RINT value") +
    my_theme
  patient.plots[[length(patient.plots) + 1]] <- p
}

combined_patient_plot <-
    wrap_plots(patient.plots, ncol = 3) +
    plot_layout(guides = "collect") +     
    plot_annotation(
        title = "Age at diagnosis-CBC correlations in T-LGLL",
        theme = theme(
        plot.title = element_text(hjust = 0.5,face = "bold",size = 8))
    )

ggsave(file.path(fig.path, "FigureS1.png"),combined_patient_plot, width = 16, height = 17, units = "cm", dpi = 300)


# ----------------- Figure S2. associations in NHANES    --------------------

nhanes.plots <- list()
normals_res  <- tibble()

for (var in vars) {
  total  <- normal |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
  total2 <- rm.outliers(total)
  total2[[var]] <- RINT(total2[[var]])

  model   <- lm(as.formula(paste0(var, " ~ Age")), data = total2)
  total.p <- coef(summary(model))[2, 4]
  normals_res <- bind_rows(normals_res, tibble(name = var, pval = total.p))

  p <- ggscatter(total2, x = "Age", y = var, add = "reg.line", conf.int = TRUE,
                  cor.coeff.args = list(method = "pearson"), alpha = 0.01) +
    coord_cartesian(xlim = c(20, 80)) +
    stat_cor(label.x = 30) +
    ggtitle(var) +
    labs(x = "Age", y = "RINT value") +
    my_theme
  nhanes.plots[[length(nhanes.plots) + 1]] <- p
}

combined_nhanes_plot <-
    wrap_plots(nhanes.plots, ncol = 3) +
    plot_layout(guides = "collect") +     
    plot_annotation(
        title = "Age at diagnosis-CBC correlations in NHANES cohort",
        theme = theme(
        plot.title = element_text(hjust = 0.5,face = "bold",size = 8))
    )

ggsave(file.path(fig.path, "FigureS2.png"),combined_nhanes_plot, width = 16, height = 17, units = "cm", dpi = 300)


# ----------------- Figure 1A-B. Cocor comparison (NHANES vs T-LGLL) --------------------

# Figure 1A. 

rows.forest  <- list()
plots <- list() 
ks.p <- tibble()
sig_cocor_res = NULL

for (var in vars) {
  normals  <- normal  |> dplyr::select(all_of(var), Age) |> mutate(Group = "Healthy") |>
    na.omit() |> as.data.frame()
  patients <- patient |> dplyr::select(all_of(var), Age) |> mutate(Group = "LGLL") |>
    na.omit() |> as.data.frame()

  normals2  <- rm.outliers(normals)
  patients2 <- patients[point.keepers(patients), ]
  
  if (nrow(patients2) < 15) next

  res.ks    <- ks.test(normals2 |> pull(var), patients2 |> pull(var))
  ks.p      <- bind_rows(ks.p, tibble(name = var, pval = res.ks$p.value))

  normals2[[var]]  <- RINT(normals2[[var]])
  patients2[[var]] <- RINT(patients2[[var]])

  res.cocor <- cocor(as.formula(paste0("~ Age + ", var, "| Age + ", var)),
                     data = list(normals2, patients2))
  ct.norm   <- cor.test(normals2[[var]],  normals2$Age,  method = "pearson")
  ct.pat    <- cor.test(patients2[[var]], patients2$Age, method = "pearson")
  pat.p     <- coef(summary(lm(as.formula(paste0(var, " ~ Age")), data = patients2)))[2, 4]

  rows.forest[[length(rows.forest) + 1]] <- bind_rows(
    tibble(var = var, group = "Healthy", r = as.numeric(ct.norm$estimate),
           lower = ct.norm$conf.int[1], upper = ct.norm$conf.int[2],
           pval = ct.norm$p.value, cocor.p = res.cocor@fisher1925$p.value,
           n = nrow(normals2)),
    tibble(var = var, group = "LGLL",    r = as.numeric(ct.pat$estimate),
           lower = ct.pat$conf.int[1],  upper = ct.pat$conf.int[2],
           pval = ct.pat$p.value,  cocor.p = res.cocor@fisher1925$p.value,
           n = nrow(patients2))
  )
  
  # get significant results only for both cocor and assoc in T-LGLL 
  if (res.cocor@fisher1925$p.value >= 0.05 || pat.p >= 0.05) next

  cat("Significant cocor result for", var, ":", res.cocor@fisher1925$p.value, "\n")
  combined.df <- bind_rows(normals2, patients2) |>
    mutate(Group = factor(Group, levels = c("Healthy", "LGLL")))

  p <- ggscatter(combined.df, x = "Age", y = var, col = "Group",
                  add = "reg.line", conf.int = TRUE,
                  cor.coeff.args = list(method = "pearson")) +
    coord_cartesian(xlim = c(20, 80)) +
    stat_cor(aes(col = Group), label.x = 30) +
    ggtitle(paste0(var, " \n(cocor p = ", format.pval(res.cocor@fisher1925$p.value, 2), ")")) +
    ylab("RINT values") + 
    my_theme

    p$layers[[1]]$aes_params <- list()
    p$layers[[1]]$mapping    <- aes(Age, .data[[var]], color = Group, alpha = Group)
    p <- p + scale_alpha_manual(values = c(0.04, 0.3), name = "Group") + 
            scale_color_manual(values = c("Healthy" = "grey55", "LGLL" = "black"), name = "Group") +
            scale_fill_manual( values = c("Healthy" = "grey55", "LGLL" = "grey30"),  name = "Group") 
            

    plots[[length(plots)+1]] <- p

}

fig1A <-
    wrap_plots(plots, ncol = 3) & 
    theme(legend.position ="none")
    
ggsave(file.path(fig.path, "Figure1A.png"),fig1A, width = 16, height = 5, units = "cm", dpi = 300)

# Figure 1B. 

forest.data <- bind_rows(rows.forest)
cocor.p.adj <- forest.data |>
  filter(group == "Healthy") |>
  mutate(cocor.p.adj = p.adjust(cocor.p, method = "BH")) |>  # BH-corrected cocor p-values
  dplyr::select(var, cocor.p.adj)

forest.data <- left_join(forest.data, cocor.p.adj, by = "var")

cocor.sig.vars <- forest.data |>
  filter(group == "Healthy", cocor.p.adj < 0.05) |>
  pull(var)

var.order <- forest.data |> filter(group == "Healthy") |> arrange(r) |> pull(var)

forest.data <- forest.data |>
  mutate(
    sig   = ifelse(pval < 0.05, "p < 0.05", "n.s."),
    var   = factor(var, levels = var.order),
    group = factor(group, levels = c("Healthy", "LGLL"))
  )

ylabels <- setNames(
  ifelse(levels(forest.data$var) %in% cocor.sig.vars,
         paste0(levels(forest.data$var), "  *"),
         levels(forest.data$var)),
  levels(forest.data$var)
)

p.forest <- ggplot(forest.data, aes(x = r, y = var, color = group, group = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#c3c2b7", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                    height = 0.25, linewidth = 0.4, position = position_dodge(0.6)) +
    geom_point(aes(shape = sig), size = 3, position = position_dodge(0.6)) +
#   scale_color_manual(values = c("Healthy" = "#2a78d6", "LGLL" = "#eb6834"), name = "Group") +
    scale_color_manual(values = c("Healthy" = "grey55", "LGLL" = "black"), name = "Group") +
    scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1), name = "Significance") +
    scale_y_discrete(labels = ylabels) +
    labs(x = "Pearson r with Age", y = NULL,
        caption = "* correlation significantly different between Healthy and LGLL (cocor BH-adjusted p < 0.05)") +
    my_theme + 
    theme_classic() + 
    guides(shape = guide_legend(order = 1,override.aes = list(size = 2)),
        color = guide_legend(order = 2,override.aes = list(size = 2,linewidth = 0.4))) +
    theme(
        axis.title.x = element_text(size= 7), 
        legend.title = element_text(
        size = 6,
        face = "plain"
        ),
        legend.text = element_text(size = 7),
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(5, "mm"),

        legend.spacing.y = unit(0.5, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 1),

        plot.caption = element_text(size = 7,color = "#52514e")
    )
    
  
ggsave(file.path(fig.path, "Figure1B.png"),p.forest, width = 16, height = 6, units = "cm", dpi = 300)


# ----------------- Figure S3. subsampling sensitivity within NHANES data --------------------

set.seed(42)
n_iter <- 1000L

subsamp.rows <- list()
for (var in vars) {
  nhanes.sub <- normal |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
  nhanes.sub <- rm.outliers(nhanes.sub)

  patients.sub <- patient |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
  if (nrow(patients.sub) < 15) next
  patients.sub <- patients.sub[point.keepers(patients.sub), ]
  if (nrow(patients.sub) < 15) next
  n_lgll <- nrow(patients.sub)
  if (nrow(nhanes.sub) < n_lgll) next

  patients.sub[[var]] <- RINT(patients.sub[[var]])
  lgll.p <- coef(summary(lm(as.formula(paste0(var, " ~ Age")), data = patients.sub)))[2, 4]
  lgll.r <- cor(patients.sub[[var]], patients.sub$Age)
  nhanes.r.full <- cor(RINT(nhanes.sub[[var]]), nhanes.sub$Age)

  sig_count <- 0L
  r_vals    <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    samp        <- nhanes.sub[sample(nrow(nhanes.sub), size = n_lgll, replace = FALSE), ]
    samp[[var]] <- RINT(samp[[var]])
    p_i <- tryCatch(
      coef(summary(lm(as.formula(paste0(var, " ~ Age")), data = samp)))[2, 4],
      error = function(e) NA_real_
    )
    r_vals[i] <- cor(samp[[var]], samp$Age)
    if (!is.na(p_i) && p_i < 0.05) sig_count <- sig_count + 1L
  }

  subsamp.rows[[length(subsamp.rows) + 1]] <- tibble(
    var                   = var,
    nhanes_n              = nrow(nhanes.sub),
    lgll_n                = n_lgll,
    nhanes_r              = nhanes.r.full,
    nhanes_r_matched_mean = mean(r_vals, na.rm = TRUE),
    nhanes_r_matched_sd   = sd(r_vals, na.rm = TRUE),
    lgll_r                = lgll.r,
    lgll_p                = lgll.p,
    emp_power             = sig_count / n_iter
  )
}

subsamp.df <- bind_rows(subsamp.rows)

  subsamp.df <- subsamp.df |>
    mutate(
      lgll_p_adj = p.adjust(lgll_p, method = "BH"),
      lgll_sig   = lgll_p_adj < 0.06,
      lgll_label = ifelse(lgll_sig, "LGLL adj.p < 0.05", "LGLL n.s.")
    )

  var.order.sp <- subsamp.df |> arrange(emp_power) |> pull(var)
  subsamp.df   <- subsamp.df |> mutate(var = factor(var, levels = var.order.sp))

  figS3 <- ggplot(subsamp.df, aes(x = emp_power, y = var)) +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "#888888", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = emp_power, y = var, yend = var),
                 color = "grey80", linewidth = 0.4) +
    geom_point(aes(color = lgll_label, shape = lgll_label), size = 3) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
    scale_color_manual(
      values = c("LGLL adj.p < 0.05" = "#0072B2", "LGLL n.s." = "#D55E00"), name = NULL
    ) +
    scale_shape_manual(
      values = c("LGLL adj.p < 0.05" = 16, "LGLL n.s." = 1), name = NULL
    ) +
    labs(
      x       = "Empirical detection rate in NHANES subsamples at matched N (1,000 draws, p < 0.05)",
      y       = NULL,
      caption = "Dashed line = 80% power threshold"
    ) +
    my_theme + 
    scale_color_manual(
      values = c("LGLL adj.p < 0.05" = "black", "LGLL n.s." = "black"), name = NULL
    ) +
    scale_shape_manual(
      values = c("LGLL adj.p < 0.05" = 16, "LGLL n.s." = 1), name = NULL
    ) +
    theme(axis.title.x = element_text(size =7),
        legend.position = "right", 
        legend.text = element_text(size = 7),
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(5, "mm"),

        legend.spacing.y = unit(0.5, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 1),

        plot.caption = element_text(size = 6, color = "grey55")) 
    

    ggsave(file.path(fig.path, "FigureS3.png"),figS3, width = 12, height = 6, units = "cm", dpi = 300)


# ----------------- Figure S4. Sex-stratified correlation forest  --------------------


sex.forest.rows <- list()
for (sx in c("M", "F")) {
  normal.sx  <- normal      |> filter(Sex == sx)
  patient.sx <- patient.sex |> filter(Sex == sx)

  for (var in vars) {
    normals2  <- normal.sx |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
    normals2  <- rm.outliers(normals2)
    normals2[[var]] <- RINT(normals2[[var]])
    ct.norm   <- cor.test(normals2[[var]], normals2$Age, method = "pearson")

    patients2 <- patient.sx |> dplyr::select(all_of(var), Age) |> na.omit() |> as.data.frame()
    if (nrow(patients2) < 10) next
    patients2 <- patients2[point.keepers(patients2), ]
    if (nrow(patients2) < 10) next
    patients2[[var]] <- RINT(patients2[[var]])
    ct.pat <- cor.test(patients2[[var]], patients2$Age, method = "pearson")

    cocor.p.sx <- cocor(as.formula(paste0("~ Age + ", var, "| Age + ", var)),
                        data = list(normals2, patients2))@fisher1925$p.value

    sex.forest.rows[[length(sex.forest.rows) + 1]] <- bind_rows(
      tibble(var = var, Sex = sx, group = "Healthy", r = as.numeric(ct.norm$estimate),
             lower = ct.norm$conf.int[1], upper = ct.norm$conf.int[2],
             pval = ct.norm$p.value, cocor.p = cocor.p.sx, n = nrow(normals2)),
      tibble(var = var, Sex = sx, group = "LGLL",    r = as.numeric(ct.pat$estimate),
             lower = ct.pat$conf.int[1],  upper = ct.pat$conf.int[2],
             pval = ct.pat$p.value,  cocor.p = cocor.p.sx, n = nrow(patients2))
    )
  }
}
sex.forest.df <- bind_rows(sex.forest.rows)

cocor.adj.sex <- sex.forest.df |>
  filter(group == "Healthy") |>
  group_by(Sex) |>
  mutate(cocor.p.adj = p.adjust(cocor.p, method = "BH")) |>
  ungroup() |>
  dplyr::select(var, Sex, cocor.p.adj)
sex.forest.df <- left_join(sex.forest.df, cocor.adj.sex, by = c("var", "Sex"))

cocor.sig.sex <- sex.forest.df |>
  filter(group == "Healthy", cocor.p.adj < 0.05) |>
  mutate(var_sex = paste0(var, "_", Sex)) |>
  pull(var_sex)

sex.forest.df <- sex.forest.df |>
  mutate(
    sig      = ifelse(pval < 0.05, "p < 0.05", "n.s."),
    var      = factor(var, levels = intersect(var.order, unique(var))),
    group    = factor(group, levels = c("Healthy", "LGLL")),
    var_sex  = paste0(var, "_", Sex),
    label_star = ifelse(var_sex %in% cocor.sig.sex, "*", "")
  )

x.max.sf  <- max(sex.forest.df$upper, na.rm = TRUE)
star.annot <- sex.forest.df |>
  filter(group == "Healthy", label_star == "*") |>
  mutate(x_pos = x.max.sf * 1.08)

p.sex.forest <- ggplot(sex.forest.df, aes(x = r, y = var, color = group, group = group)) +
  facet_wrap(~Sex, scales = "fixed",
             labeller = labeller(Sex = c("M" = "Male", "F" = "Female"))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#c3c2b7", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.25, linewidth = 0.6, position = position_dodge(0.6)) +
  geom_point(aes(shape = sig), size = 3, position = position_dodge(0.6)) +
  geom_text(data = star.annot, aes(x = x_pos, y = var, label = "*"),
            color = "#0b0b0b", size = 4, hjust = 0, inherit.aes = FALSE) +
  scale_color_manual(values = c("Healthy" = "#2a78d6", "LGLL" = "#eb6834"), name = "Group") +
  scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1), name = "Significance") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Pearson r with Age", y = NULL,
       caption = "* cocor p < 0.05 (BH-corrected): correlation differs between Healthy and LGLL") +
  scale_color_manual(values = c("Healthy" = "grey55", "LGLL" = "black"), name = "Group") +
  my_theme + 
  theme(axis.title.x = element_text(size= 7),
        # axis.text.y     = element_text(size = 7),
        # strip.text      = element_text(size = 8, face = "bold"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(5, "mm"),

        legend.spacing.y = unit(0.5, "mm"),
        legend.box.spacing = unit(1, "mm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 1),
        plot.caption    = element_text(size = 7, color = "#52514e"))  


ggsave(file.path(fig.path, "FigureS4.png"),p.sex.forest, width = 16, height = 6, units = "cm", dpi = 300)


# ----------------- Figure 1C. Age-adjusted z-score (vs NHANES) --------------------

patient.sex <- patient |>
  mutate(Sex = cli.new$Sex)

zscore.rows <- list()
for (var in vars) {
  nhanes.sub <- normal |> dplyr::select(all_of(var), Age, Sex) |> na.omit() |> as.data.frame()
  nhanes.sub <- rm.outliers(nhanes.sub)

  model.z  <- lm(as.formula(paste0(var, " ~ Age + Sex")), data = nhanes.sub)
  sigma.z  <- sigma(model.z)

  pat.sub <- patient.sex |> dplyr::select(all_of(var), Age, Sex, STAT3) |> na.omit() |>
    as.data.frame()
  if (nrow(pat.sub) < 10) next

  predicted.z <- predict(model.z, newdata = pat.sub)
  z           <- (pat.sub[[var]] - predicted.z) / sigma.z

  zscore.rows[[length(zscore.rows) + 1]] <- tibble(
    var   = var,
    z     = z,
    Sex   = pat.sub$Sex,
    STAT3 = pat.sub$STAT3,
    Age = pat.sub$Age
  )
}

zscore.df <- bind_rows(zscore.rows)

var.order.z <- zscore.df |>
  group_by(var) |>
  summarise(med = median(z, na.rm = TRUE), .groups = "drop") |>
  arrange(med) |> pull(var)

zscore.df <- zscore.df |> mutate(var = factor(var, levels = var.order.z))

wilcox.res <- zscore.df |>
  group_by(var) |>
  summarise(p.raw = wilcox.test(z, mu = 0)$p.value, .groups = "drop") |>
  mutate(p.adj = p.adjust(p.raw, method = "BH"),
         sig   = ifelse(p.adj < 0.05, "p < 0.05", "n.s."))

zscore.df <- left_join(zscore.df, wilcox.res |> dplyr::select(var, sig), by = "var")

# p.zscore <- ggplot(zscore.df, aes(x = var, y = z)) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "#c3c2b7", linewidth = 0.5) +
#     geom_boxplot(aes(fill = sig), outlier.shape = NA, alpha = 0.5, width = 0.5) +
#     geom_jitter(width = 0.2, alpha = 0.25, size = 0.5, color = "#52514e") +
#     scale_fill_manual(values = c("p < 0.05" = "grey55", "n.s." = "grey75"),
#                         name = "Wilcoxon vs 0\n(BH-corrected)") +
#     labs(x = NULL, y = "Age- and sex-adjusted z-score\n(T-LGLL patient vs NHANES reference)") +
#     coord_flip(ylim = c(-10, 10)) + 
#     my_theme + 
#     theme_classic() +
#     theme(legend.position = "none", 
#         axis.title.x = element_text(size =7)) 

    p.zscore <- ggplot(zscore.df,aes(x = z, y = var)) +
        geom_vline(xintercept = 0,linetype = "dashed",color = "#c3c2b7",linewidth = 0.5) +
        geom_boxplot(fill = "grey55",orientation = "y",outlier.shape = NA,alpha = 0.5,width = 0.5) +
        geom_jitter(width = 0,height = 0.2,alpha = 0.25,size = 0.5,color = "#52514e") +
        scale_y_discrete(limits = var.order.z,drop = FALSE) +
        coord_cartesian(xlim = c(-10, 10)) +
        labs(x = paste0("Age- and sex-adjusted z-score\n","(T-LGLL patient vs NHANES reference)"),
            y = NULL) +
        my_theme +
        theme_classic() +
        theme(
        legend.position = "none",
        axis.title.x = element_text(size = 7)
    )

    ggsave(file.path(fig.path, "Figure1C.png"),p.zscore, width = 16, height = 6, units = "cm", dpi = 300)


# STAT3 association 

zscore.stat3 <- zscore.df |>
  mutate(
    STAT3 = factor(STAT3),
    STAT3 = relevel(STAT3, ref = "WT")
  ) |>
  filter(
    !is.na(z),
    !is.na(Age),
    !is.na(STAT3)
  )
model_add <- lm(z ~ Age + STAT3, data = sub)
model_int <- lm(z ~ Age * STAT3, data = sub)

anova(model_add, model_int)


# ----------------- section 2. time to treatment related. but see whether the CBC parameters can predict the treatment need  --------------------

#  regression: CBC parameters predict treatment need

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