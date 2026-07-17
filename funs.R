my_theme <- theme_minimal(base_size = 8) +
  theme(
    plot.title      = element_text(size = 10),
    axis.title      = element_text(size = 10),
    axis.text       = element_text(size = 8),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8),
    strip.text      = element_text(size = 9),
    text            = element_text(family = "Arial")
  )

clinical.path <- file.path("data", "clinical.csv")
germline.variants.path <- file.path("data", "germline_variant_summary.tsv.gz")
germline.patients.path <- file.path("data", "germline_patient_summary.tsv.gz")
somatic.path <- file.path("data", "somatic_mutations.tsv.gz")
rnaseq.path <- file.path("data", "counts.tsv.gz")

dice.metadata.path <- file.path("data", "DICE", "mmc1.xlsx")
dice.cd8.path <- file.path("data", "DICE", "CD8_NAIVE_TPM.csv")

nhanes.cbc.path <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/CBC_L.xpt"
nhanes.demo.path <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.xpt"
nhanes.cols.path <- file.path("data", "CDC_corresponding_columns.xlsx")

single_column_width <- 8
onepointfive_column_width <- 12
two_column_width <- 16
