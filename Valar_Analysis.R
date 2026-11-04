set.seed(111)

# ================================= LIBRARIES ==================================

library(tidyverse)
library(openxlsx)
library(edgeR)
library(proDA)
library(fgsea)
library(msigdbr)
library(ggh4x)
library(patchwork)
library(ragg)
library(scales)
library(png)
library(grid)

# ================================ PARAMETERS ==================================

input_file_RNA_counts <-
  "./Data/RNA/salmon.merged.gene_counts.tsv"

input_file_RNA_metadata <-
  "./Data/RNA/BIFSAMPLE.xlsx"

input_file_LV <-
  "./Data/Proteomics/P03_directLFQ_unfiltered_LV.tsv"

input_file_HT <-
  "./Data/Proteomics/P03_directLFQ_unfiltered_HT.tsv"

input_file_SampleInfo <-
  "./Data/Proteomics/P03_metadata_sdrf.tsv"

input_file_GeneSymbolMapping <-
  c(LV = "./Data/Proteomics/idmapping_2026_04_27_liver.tsv",
    HT = "./Data/Proteomics/idmapping_2026_04_27_heart.tsv")

# fgsea results used in the manuscript (fgsea p values and NES are permutation
# based, so these frozen tables are the inputs of record for figures and tables)

input_file_enrichment <-
  c(Br = "./Data/Enrichment/cortex_ALL.xlsx",
    Hr = "./Data/Enrichment/heart_ALL.xlsx",
    Lv = "./Data/Enrichment/liver_ALL.xlsx",
    LV_protein = "./Data/Enrichment/MSIG_REACTOME_Lv_No_4W.xlsx")

res_folder <-
  "./Results/"

fig_folder <-
  "./Results/Figures/"

tab_folder <-
  "./Results/Tables/"

# The frozen tables were
# computed with msigdbr 10.0.2 + msigdbdf 24.1.1; msigdbr 26.1.1 ships different gene
# sets (1,221 instead of 1,143 sets of 15-500 genes), so a faithful re-run needs
# those versions.

recompute_enrichment = FALSE

# Proteomics: tissue-specific mean CV threshold between technical replicates
cv_threshold = c(LV = 0.5, HT = 0.4)

dir.create(fig_folder, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_folder, showWarnings = FALSE, recursive = TRUE)

options(width = 200)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                    ---- RNA SEQUENCING (limma-voom) ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ---- Reading data ----

T = read.delim(input_file_RNA_counts)

M = read.xlsx(input_file_RNA_metadata)

# ---- Cleaning and unifying data sets ----

gene_names = T %>%
  select(gene_id, gene_name) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id))

meta = M %>%
  select(sample = `Sample.Name.in.Report.(Required)`,
         Group.Name2) %>%
  distinct()

colnames(meta)[2] = "region"

T_long = T %>%
  pivot_longer(cols = -c(gene_id, gene_name),
               names_to = "sample",
               values_to = "n") %>%
  left_join(meta, by = "sample") %>%
  mutate(region = recode(region,
                         "Liver"  = "Lv",
                         "Heart"  = "Hr",
                         "Cortex" = "Br",
                         "Lung"   = "Lu")) %>%
  extract(col = sample,
          into = c("condition", "replicate"),
          regex = "^(.*)_(\\d+)$",
          remove = FALSE) %>%
  mutate(replicate = as.integer(replicate)) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id))

# ---- Subsetting dataset ----

T_long = T_long %>%
  filter(condition != "Sh_coal",
         condition != "Sh_el",
         replicate < 27) %>%
  select(-sample) %>%
  arrange(region, condition, replicate, gene_id)

# ---- Pivoting dataset ----

T_spread = T_long %>%
  select(condition, replicate, region, gene_id, n) %>%
  pivot_wider(names_from = gene_id,
              values_from = n)

T_spread = as.data.frame(T_spread)

metadata_idxs = 1:3

# ---- DEGS with limma voom ----

DEGS = data.frame()

for (reg in c("Br", "Lv", "Hr")) {

  cat("\n")
  cat(reg)
  cat("\n")

  D = T_spread %>% filter(region == reg)

  cd = t(D[, 4:ncol(D)])
  colnames(cd) = paste(D$region, D$condition, D$replicate, sep = "_")

  md = D[, metadata_idxs]
  rownames(md) = paste(D$region, D$condition, D$replicate, sep = "_")

  md$condition = factor(x = md$condition,
                        levels = c("Ctr", "Noise"))

  dge = DGEList(counts = cd)

  dge$samples = cbind(dge$samples, md)

  design = model.matrix(~ condition, data = dge$samples)

  keep = filterByExpr(dge, design)
  dge = dge[keep, , keep.lib.sizes = FALSE]

  dge = calcNormFactors(dge)

  # adaptive.span = FALSE: voom default of limma < 3.68 used for the manuscript
  v = voom(dge, design, adaptive.span = FALSE)
  fit = lmFit(v, design)
  fit = eBayes(fit)

  contrast = makeContrasts(contrasts = "conditionNoise", levels = design)
  fit2 = contrasts.fit(fit, contrast)
  fit2 = eBayes(fit2)
  res = topTable(fit2, n = Inf, sort.by = "none")

  res$Region = reg
  res$ident1 = "Noise"
  res$ident2 = "Ctr"
  res$gene_id = rownames(res)

  DEGS = rbind(DEGS, res)

}

# ---- Add gene symbols to DEGs ----

DEGS$gene_symbol = gene_names$gene_name[match(DEGS$gene_id, gene_names$gene_id)]

rownames(DEGS) = NULL

write.xlsx(DEGS, file = paste0(res_folder, "RNA_DEGS_unfiltered.xlsx"))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                 ---- PROTEOMICS: LIVER AND HEART (proDA) ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ---- Importing data-sets ----

P_LV = read.csv(input_file_LV, sep = "\t") %>% pivot_longer(cols = 2:41,
                                                          names_to = "Sample",
                                                          values_to = "n")
colnames(P_LV)[1] = "protein"
P_LV$Region = "LV"

P_HT = read.csv(input_file_HT, sep = "\t") %>% pivot_longer(cols = 2:41,
                                                          names_to = "Sample",
                                                          values_to = "n")
colnames(P_HT)[1] = "protein"
P_HT$Region = "HT"

P_raw = rbind(P_LV, P_HT)

P_raw = P_raw %>%
  separate(Sample, into = c("Sample", "Replicate"), sep = "_") %>%
  mutate(Replicate = as.numeric(Replicate)) %>%
  mutate(Sample = paste0(Sample, "_", Replicate)) %>%
  select(-Replicate)

# ---- Assigning metadata ----

SampleInfo = read.csv(input_file_SampleInfo, sep = "\t")
SampleInfo$Region = recode(SampleInfo$characteristics.organism.part.,
                           "liver" = "LV",
                           "heart" = "HT")
SampleInfo$SampleReplicate = paste(SampleInfo$source.name,
                                   SampleInfo$comment.technical.replicate.,
                                   sep = "_")

P_raw$Condition = SampleInfo$characteristics.disease.[match(paste(P_raw$Region, P_raw$Sample),
                                                            paste(SampleInfo$Region, SampleInfo$SampleReplicate))]
P_raw = P_raw %>% separate(Sample,
                           into = c("Sample", "Replicate"),
                           sep = "_")
P_raw$n[P_raw$n == 0] = NA

P_raw$Condition = str_replace(P_raw$Condition, pattern = "AN_4 days", replacement = "4D")
P_raw$Condition = str_replace(P_raw$Condition, pattern = "AN_28 days", replacement = "4W")
P_raw$Condition = str_replace(P_raw$Condition, pattern = "CTRL", replacement = "Control")

# ---- Average CoV by protein (replicate level) ----

P_CoV = P_raw %>%
  group_by(Region, protein, Sample, Condition) %>%
  summarise(
    cv_n = ifelse(sum(!is.na(n)) > 1, sd(n, na.rm = TRUE) / mean(n, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  group_by(Region, protein) %>%
  summarise(
    mean_cv = mean(cv_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cv_threshold = cv_threshold[Region])

p = ggplot(P_CoV, aes(x = mean_cv)) +
  geom_histogram(bins = 100) +
  geom_vline(aes(xintercept = cv_threshold), color = "red", linewidth = 1, linetype = "dashed") +
  facet_wrap(~ Region, scales = "free_y") +
  labs(x = "Mean CoV between technical replicates", y = "Count") +
  theme_minimal()

ggsave(paste0(res_folder, "QC_Proteomics_replicate_CoV.pdf"), p, width = 8, height = 4)

# ---- Filter proteins that exhibit high measurement variation between technical replicates ----

to_remove = P_CoV %>%
  filter(mean_cv > cv_threshold) %>%
  select(Region, protein)

P_raw = P_raw %>% anti_join(to_remove, by = c("Region", "protein"))

print(table(to_remove$Region))

# ---- Calculate mean quantification between technical replicates ----
#(in the case of one numeric value and one NA value, keep the numeric value)

P = P_raw %>%
  group_by(Region, protein, Sample, Condition) %>%
  summarise(n = if (all(is.na(n))) {
    NA_real_
  } else if (any(is.na(n))) {
    n[!is.na(n)][1]  # keep the non-NA value
  } else {
    mean(n)  # take mean of two values
  }, .groups = "drop")

# ---- Removing entries that exhibit too many NAs ----
#(within one tissue, for each protein, if we observe that in all conditions
#there are >=30% missing values, we remove those quantification)

P_missing = P %>%
  group_by(Region, protein, Condition) %>%
  summarise(na_prop = mean(is.na(n)), .groups = "drop") %>%
  group_by(Region, protein) %>%
  summarise(remove = all(na_prop >= 0.3), .groups = "drop")

print(table(P_missing$Region, P_missing$remove))

P = P %>%
  left_join(P_missing, by = c("Region", "protein")) %>%
  filter(!remove) %>%
  select(-remove)

# ---- DAPS with proDA ----

DAPS = data.frame()

for (region in c("LV", "HT")) {

  cat("\n")
  cat(region)
  cat("\n")

  P_spread = P %>%
    filter(Region == region) %>%
    select(-Region) %>%
    pivot_wider(names_from = protein, values_from = n)

  expr_matrix = P_spread %>%
    select(-Condition, -Sample) %>%
    as.matrix() %>%
    t() %>%
    log2()

  expr_matrix = expr_matrix[which(rowSums(is.na(expr_matrix)) != dim(expr_matrix)[2]), ]

  expr_matrix = proDA::median_normalization(expr_matrix)

  sample_info = P_spread %>%
    select(Condition)

  sample_info$Condition = factor(sample_info$Condition,
                                 levels = c("Control", "4D", "4W"))

  fit = proDA(
    data = expr_matrix,
    design = ~ Condition,
    col_data = sample_info,
    reference_level = "Control"
  )

  DEG_4D = test_diff(fit, "Condition4D")
  DEG_4D$ident.1 = "4D"
  DEG_4D$ident.2 = "Control"

  DEG_4W = test_diff(fit, "Condition4W")
  DEG_4W$ident.1 = "4W"
  DEG_4W$ident.2 = "Control"

  temp = rbind(DEG_4D, DEG_4W)

  # Gene names mapping (unmapped accessions are kept as they are)

  ProteinToGene = read.csv(input_file_GeneSymbolMapping[region], sep = "\t")
  protein_to_gene = setNames(ProteinToGene$To, ProteinToGene$From)

  gene_map = temp %>%
    distinct(name) %>%
    mutate(protein_id = name) %>%
    separate_rows(protein_id, sep = ";") %>%
    mutate(gene = unname(protein_to_gene[protein_id]),
           gene = ifelse(is.na(gene), protein_id, gene)) %>%
    group_by(name) %>%
    summarise(GeneSymbol = paste(gene, collapse = ";"), .groups = "drop")

  temp = temp %>%
    left_join(gene_map, by = "name") %>%
    arrange(adj_pval)

  temp$Region = region

  DAPS = rbind(DAPS, temp)

}

DAPS_LV = DAPS %>% filter(Region == "LV") %>% select(-Region)
DAPS_HT = DAPS %>% filter(Region == "HT") %>% select(-Region)

write.xlsx(DAPS_LV, file = paste0(res_folder, "Proteomics_DiffExpr_ProDA_Liver.xlsx"))
write.xlsx(DAPS_HT, file = paste0(res_folder, "Proteomics_DiffExpr_ProDA_Heart.xlsx"))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#              ---- GENE-SET ENRICHMENT (fgsea, Hallmark + Reactome) ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ---- Gene sets (MSigDB Hallmark + Reactome, mouse, 15-500 genes) ----

msig_h = msigdbr(species = "Mus musculus", collection = "H") %>%
  select(gs_name, gene_symbol)

msig_reactome = msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME") %>%
  select(gs_name, gene_symbol)

msig = bind_rows(msig_h, msig_reactome)

msig_sets = msig %>%
  group_by(gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  mutate(size = lengths(genes)) %>%
  filter(size >= 15, size <= 500)

msig_all = setNames(msig_sets$genes, msig_sets$gs_name)

cat("\nGene sets with 15-500 genes:", length(msig_all), "\n")

# ---- Optional: recompute enrichment ----
#(RNA: genes ranked by the moderated t; protein: proteins ranked by the proDA t
#statistic of 4W vs Control; duplicated symbols keep the largest |t|)

if (recompute_enrichment) {

  dir.create(paste0(res_folder, "Enrichment_recomputed"), showWarnings = FALSE)

  for (reg in c("Br", "Hr", "Lv", "LV_protein")) {

    if (reg == "LV_protein") {
      ranks_tbl = DAPS_LV %>%
        filter(ident.1 == "4W", !is.na(GeneSymbol), !is.na(t_statistic)) %>%
        transmute(gene_symbol = GeneSymbol, t = t_statistic)
    } else {
      ranks_tbl = DEGS %>%
        filter(Region == reg, !is.na(gene_symbol), !is.na(t)) %>%
        select(gene_symbol, t)
    }

    ranks_tbl = ranks_tbl %>%
      group_by(gene_symbol) %>%
      slice_max(abs(t), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(t))

    ranks = setNames(ranks_tbl$t, ranks_tbl$gene_symbol)

    res = fgsea(pathways = msig_all, stats = ranks) %>%
      as_tibble() %>%
      arrange(padj, desc(abs(NES)))

    write.xlsx(res, paste0(res_folder, "Enrichment_recomputed/", reg, "_ALL.xlsx"))

    frozen = read.xlsx(input_file_enrichment[reg]) %>%
      inner_join(res %>% select(pathway, NES_new = NES, padj_new = padj), by = "pathway")

    cat(sprintf("%s: %d gene sets; max |dNES| vs manuscript %.3f; FDR<0.05 agreement %d / %d\n",
                reg, nrow(frozen), max(abs(frozen$NES - frozen$NES_new)),
                sum((frozen$padj < 0.05) == (frozen$padj_new < 0.05)), nrow(frozen)))
  }
}

# ---- Load enrichment results used in the manuscript ----

E_cortex = read.xlsx(input_file_enrichment["Br"])
E_heart  = read.xlsx(input_file_enrichment["Hr"])
E_liver  = read.xlsx(input_file_enrichment["Lv"])
P_liver  = read.xlsx(input_file_enrichment["LV_protein"])

stopifnot(nrow(E_cortex) == 1143, nrow(E_heart) == 1143, nrow(E_liver) == 1143)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                       ---- MANUSCRIPT PLOTS ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ---- Shared style ----

fig_width_mm = 180

figX_height_A_mm = 108

figX_height_Bcurve_mm    = 30  # running enrichment score
figX_height_Bticks_mm    = 12  # barcode + rank axis
figX_height_Blollipop_mm = 50  # fold-change plot, incl. its gene labels
figX_height_BC_mm = figX_height_Bcurve_mm + figX_height_Bticks_mm +
  figX_height_Blollipop_mm - 5

figX_height_mm = figX_height_A_mm + figX_height_BC_mm

tilt_angle = 60

figX_width_B_mm = 0.60 * fig_width_mm
figX_width_C_mm = fig_width_mm - figX_width_B_mm

fig7_height_A_mm = 62

col_nes_low  = "#008585"
col_nes_mid  = "#fbf2c4"
col_nes_high = "#c7522a"

col_direction = c(
  "Down in both" = col_nes_low,
  "Up in both"   = col_nes_high,
  "Opposite"     = "#6A51A3"
)

nes_limits = c(-3.5, 3.5)
nlp_cap    = 15
nlp_breaks = c(2, 5, 10, 15)
size_range = c(1.2, 5.2)

fdr_levels = c("FDR < 0.05", "FDR ≥ 0.05")

theme_panel = theme_classic(base_size = 9) +
  theme(
    plot.caption = element_text(size = 6.5, colour = "grey35", hjust = 0),
    plot.caption.position = "plot",
    plot.margin = margin(t = 3, r = 4, b = 2, l = 3, unit = "mm"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.5, face = "bold"),
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    panel.spacing.y = unit(1.5, "mm"),
    legend.key.size = unit(3.5, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

# ---- FIGURE 6A: cortex, heart and liver RNA pathway enrichment ----
#(colour = NES, size = -log10(FDR) capped at 15, filled = FDR < 0.05;
#ref_nes / ref_fdr are the values quoted in the manuscript, checked below)

A_spec = tribble(
  ~tissue, ~category, ~pathway, ~label, ~ref_nes, ~ref_fdr,
  "Cortex", "Inflammation", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "TNFα signaling via NF-κB (H)", 2.88, 1.1e-18,
  "Cortex", "Inflammation", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "Interferon-γ response (H)", 2.83, 2.2e-18,
  "Cortex", "Energy metabolism", "REACTOME_INTEGRATION_OF_ENERGY_METABOLISM",
  "Integration of energy metabolism (R)", -1.88, 7.2e-04,
  "Heart", "Mitochondria", "REACTOME_COMPLEX_I_BIOGENESIS",
  "Complex I biogenesis (R)", 3.10, 2.1e-21,
  "Heart", "Mitochondria", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "Oxidative phosphorylation (H)", 3.13, 1.9e-34,
  "Heart", "Mitochondria", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "TCA cycle (R)", 1.85, 6.1e-03,
  "Heart", "Fuel metabolism", "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "Mitochondrial fatty-acid β-oxidation (R)", 1.64, 0.040,
  "Heart", "Fuel metabolism", "REACTOME_GLUCOSE_METABOLISM",
  "Glucose metabolism (R)", -1.74, 3.8e-03,
  "Heart", "Fuel metabolism", "REACTOME_INTEGRATION_OF_ENERGY_METABOLISM",
  "Integration of energy metabolism (R)", -1.77, 1.2e-03,
  "Heart", "Fuel metabolism", "HALLMARK_MTORC1_SIGNALING",
  "mTORC1 signaling (H)", 1.39, 0.040,
  "Heart", "Stress response", "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES",
  "Detoxification of reactive oxygen species (R)", 2.45, 4.2e-07,
  "Heart", "Stress response", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "TNFα signaling via NF-κB (H)", 1.50, 0.018,
  "Liver", "Lipid & sterol", "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "Cholesterol biosynthesis (R)", -1.88, 8.0e-03,
  "Liver", "Lipid & sterol", "HALLMARK_FATTY_ACID_METABOLISM",
  "Fatty-acid metabolism (H)", -2.70, 6.4e-14,
  "Liver", "Lipid & sterol", "HALLMARK_PEROXISOME",
  "Peroxisome (H)", -1.93, 2.9e-04,
  "Liver", "Lipid & sterol", "HALLMARK_BILE_ACID_METABOLISM",
  "Bile-acid metabolism (H)", -1.88, 9.0e-04,
  "Liver", "Lipid & sterol", "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",
  "Activation of gene expression by SREBF/SREBP (R)", -1.10, 0.47,
  "Liver", "mTOR", "HALLMARK_MTORC1_SIGNALING",
  "mTORC1 signaling (H)", -1.93, 2.9e-05,
  "Liver", "Inflammation", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "TNFα signaling via NF-κB (H)", 2.33, 1.8e-08
)

A_category_levels = c("Lipid & sterol", "mTOR", "Inflammation", "Energy metabolism",
                      "Mitochondria", "Fuel metabolism", "Stress response")

A_dat = A_spec %>%
  left_join(bind_rows(E_cortex %>% mutate(tissue = "Cortex"),
                      E_heart %>% mutate(tissue = "Heart"),
                      E_liver %>% mutate(tissue = "Liver")) %>%
              select(tissue, pathway, NES, padj, size),
            by = c("tissue", "pathway"))

stopifnot(!any(is.na(A_dat$NES)))
stopifnot(all(abs(round(A_dat$NES, 2) - A_dat$ref_nes) < 1e-9 &
                abs(signif(A_dat$padj, 2) - A_dat$ref_fdr) < 1e-9 * pmax(A_dat$ref_fdr, 1e-30)))

A_dat = A_dat %>%
  mutate(
    nlp = -log10(padj),
    fdr = factor(ifelse(padj < 0.05, fdr_levels[1], fdr_levels[2]), levels = fdr_levels),
    tissue = factor(tissue, levels = unique(A_spec$tissue)),
    category = factor(category, levels = A_category_levels),
    row_key = paste(tissue, pathway, sep = "::"),
    row_key = factor(row_key, levels = rev(paste(A_spec$tissue, A_spec$pathway, sep = "::")))
  )

cat("\n==== Figure 6A: plotted values ====\n")
print(as.data.frame(A_dat %>% transmute(tissue, category, pathway, NES = round(NES, 4),
                                        FDR = signif(padj, 3), size, fdr)), row.names = FALSE)

A_labels = setNames(A_dat$label, as.character(A_dat$row_key))

panel_A = ggplot(A_dat, aes(x = NES, y = row_key)) +
  geom_blank() +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_point(data = filter(A_dat, padj < 0.05),
             aes(size = pmin(nlp, nlp_cap), fill = NES, shape = fdr),
             colour = "grey20", stroke = 0.3) +
  geom_point(data = filter(A_dat, padj >= 0.05),
             aes(size = pmin(nlp, nlp_cap), colour = NES, shape = fdr),
             fill = "white", stroke = 1.0) +
  scale_fill_gradient2(low = col_nes_low, mid = col_nes_mid, high = col_nes_high,
                       midpoint = 0, limits = nes_limits, oob = squish,
                       breaks = c(-3, 0, 3), name = "NES") +
  scale_colour_gradient2(low = col_nes_low, mid = col_nes_mid, high = col_nes_high,
                         midpoint = 0, limits = nes_limits, oob = squish,
                         guide = "none") +
  scale_size(range = size_range, limits = c(0, nlp_cap), breaks = nlp_breaks,
             labels = c(head(nlp_breaks, -1), paste0("≥", nlp_cap)),
             name = expression(-log[10] * "(FDR)")) +
  scale_shape_manual(values = setNames(c(21, 21), fdr_levels), drop = FALSE,
                     name = NULL) +
  guides(
    fill = guide_colourbar(order = 1, direction = "horizontal",
                           barwidth = unit(16, "mm"),
                           barheight = unit(2.5, "mm"), title.vjust = 0.9),
    size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey70",
                                                       colour = "grey20",
                                                       stroke = 0.3)),
    shape = guide_legend(order = 3, override.aes = list(
      size = 3, fill = c("grey55", "white"), colour = c("grey20", "grey40"),
      stroke = c(0.3, 1.0)))
  ) +
  scale_x_continuous(limits = nes_limits, breaks = seq(-3, 3, 1.5),
                     labels = c("-3", "", "0", "", "3")) +
  scale_y_discrete(position = "right", labels = A_labels) +
  facet_nested(tissue + category ~ ., scales = "free_y", space = "free_y",
               switch = "y",
               nest_line = element_line(colour = "grey30", linewidth = 0.4),
               strip = strip_nested(
                 text_y = elem_list_text(face = c("bold", "plain"),
                                         size = c(8, 7),
                                         angle = c(0, 0), hjust = c(1, 1)),
                 by_layer_y = TRUE, size = "variable")) +
  coord_cartesian(clip = "off") +
  labs(x = "NES", y = NULL) +
  theme_panel +
  theme(legend.position = "right",
        legend.justification = c(0, 0.5),
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(0.5, "mm"),
        legend.margin = margin(0, 2, 0, 0, unit = "mm"),
        legend.box.margin = margin(0, 0, 0, 4, unit = "mm"),
        plot.margin = margin(t = 3, r = 3, b = 2, l = 1, unit = "mm"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
        panel.spacing.y = unit(1.6, "mm"),
        axis.text.y.right = element_text(size = 6.8, hjust = 0),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave(paste0(fig_folder, "Figure6A.pdf"), panel_A,
       width = fig_width_mm / 25.4, height = figX_height_A_mm / 25.4, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(fig_folder, "Figure6A.png"), panel_A,
       width = fig_width_mm / 25.4, height = figX_height_A_mm / 25.4, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# ---- FIGURE 6B: complex I biogenesis in heart, leading edge ----
#(running enrichment score, barcode of the gene-set members with the leading
#edge in orange, and log2FC of the 41 leading-edge genes, filled = FDR < 0.05)

B_pathway = "REACTOME_COMPLEX_I_BIOGENESIS"
B_label   = "Complex I biogenesis (R)"

DEGS_Hr = DEGS %>% filter(Region == "Hr", ident1 == "Noise", ident2 == "Ctr")

# Duplicated symbols resolved by max |t|, as for the fgsea ranking
DEGS_Hr_unique = DEGS_Hr %>%
  filter(!is.na(t)) %>%
  group_by(gene_symbol) %>%
  slice_max(abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

B_res = E_heart %>% filter(pathway == B_pathway)
B_le  = str_trim(str_split(B_res$leadingEdge, ",")[[1]])

B_genes = DEGS_Hr_unique %>%
  filter(gene_symbol %in% B_le) %>%
  select(gene_symbol, logFC, adj.P.Val) %>%
  arrange(desc(logFC)) %>%
  mutate(fdr = factor(ifelse(adj.P.Val < 0.05, fdr_levels[1], fdr_levels[2]), levels = fdr_levels),
         gene_symbol = factor(gene_symbol, levels = gene_symbol))

B_sig = B_genes %>% filter(adj.P.Val < 0.05)
B_ns  = B_genes %>% filter(adj.P.Val >= 0.05)

stopifnot(nrow(B_genes) == 41, nrow(B_sig) == 29, all(B_genes$logFC > 0))
stopifnot(all.equal(round(c(range(B_genes$logFC), median(B_genes$logFC)), 2), c(0.31, 1.15, 0.67)))
stopifnot(all.equal(round(c(range(B_sig$logFC), median(B_sig$logFC)), 2), c(0.51, 1.15, 0.75)))
stopifnot(all.equal(round(range(B_ns$adj.P.Val), 3), c(0.050, 0.123)))

cat("\n==== Figure 6B: ", B_label, " ====\n", sep = "")
cat(sprintf("NES %+.2f, FDR %s; set size %d, leading edge %d\n",
            B_res$NES, formatC(B_res$padj, format = "e", digits = 1), B_res$size, length(B_le)))
print(as.data.frame(B_genes %>% transmute(gene_symbol = as.character(gene_symbol),
                                          log2FC = round(logFC, 4), FDR = signif(adj.P.Val, 3), fdr)),
      row.names = FALSE)

# Running enrichment score (as in fgsea, gseaParam = 1)

ranks_heart = sort(setNames(DEGS_Hr_unique$t, DEGS_Hr_unique$gene_symbol),
                   decreasing = TRUE)

set_genes = msig_reactome %>%
  filter(gs_name == B_pathway) %>%
  pull(gene_symbol) %>%
  unique()

stopifnot(sum(set_genes %in% names(ranks_heart)) == B_res$size)

B_curve = tibble(gene = names(ranks_heart), t = unname(ranks_heart)) %>%
  mutate(
    rank = row_number(),
    hit = gene %in% set_genes,
    step = ifelse(hit, abs(t) / sum(abs(t[hit])), -1 / sum(!hit)),
    running_es = cumsum(step)
  )

B_peak = B_curve %>% slice_max(abs(running_es), n = 1, with_ties = FALSE)
B_le_curve = B_curve %>% filter(hit, rank <= B_peak$rank) %>% pull(gene)

stopifnot(abs(B_peak$running_es - B_res$ES) < 1e-8)
stopifnot(setequal(B_le_curve, B_le))
stopifnot(abs(calcGseaStat(ranks_heart, which(names(ranks_heart) %in% set_genes)) - B_peak$running_es) < 1e-8)

B_curve = B_curve %>% mutate(in_le = gene %in% B_le)
B_hits = B_curve %>% filter(hit)

cat(sprintf("ES %.4f at rank %d of %d (recomputed; matches enrichment table)\n",
            B_peak$running_es, B_peak$rank, nrow(B_curve)))

# Bottom: leading-edge fold-changes

panel_B_lollipop = ggplot(B_genes, aes(x = gene_symbol, y = logFC)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_segment(aes(y = 0, yend = logFC, xend = gene_symbol),
               colour = "grey75", linewidth = 0.3) +
  geom_point(aes(fill = fdr), shape = 21, size = 1.5, colour = col_nes_high,
             stroke = 0.5) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, size = 2.3,
           colour = "grey15",
           label = sprintf("%d leading-edge genes", nrow(B_genes))) +
  scale_fill_manual(values = setNames(c(col_nes_high, "white"), fdr_levels),
                    drop = FALSE, name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)),
                     breaks = pretty_breaks(n = 3)) +
  labs(x = NULL, y = expression("log"[2] * " fold-change")) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  theme_panel +
  theme(axis.text.x = element_text(size = 5.3, face = "italic", angle = tilt_angle,
                                   hjust = 1, vjust = 1),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.length.x = unit(0.6, "mm"),
        legend.position = "inside",
        legend.position.inside = c(1, 0.86),
        legend.justification.inside = c(1, 1),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(2.5, "mm"),
        plot.margin = margin(t = 0, r = 3, b = 1, l = 3, unit = "mm"))

# Top: running enrichment score

panel_B_curve = ggplot(B_curve, aes(x = rank, y = running_es)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = B_peak$rank, colour = "grey55", linewidth = 0.3,
             linetype = "dashed") +
  geom_line(colour = col_nes_high, linewidth = 0.6) +
  annotate("text", x = nrow(B_curve), y = B_peak$running_es, hjust = 1, vjust = 1,
           size = 2.3, colour = "grey15",
           label = sprintf("%s, heart · NES %+.2f, FDR %s",
                           B_label, B_res$NES, formatC(B_res$padj, format = "e", digits = 1))) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = NULL, y = "Running\nenrichment score") +
  theme_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.margin = margin(t = 2, r = 3, b = 0, l = 3, unit = "mm"))

# Middle: barcode of gene-set members (carries the rank axis)

panel_B_ticks = ggplot(B_hits) +
  geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1, colour = in_le),
               linewidth = 0.25) +
  scale_colour_manual(values = c(`TRUE` = col_nes_high, `FALSE` = "grey60"),
                      labels = c(`TRUE` = "Leading edge", `FALSE` = "Other set genes"),
                      breaks = c("TRUE", "FALSE"), name = NULL) +
  scale_x_continuous(limits = range(B_curve$rank), expand = expansion(mult = 0.01),
                     breaks = seq(0, 15000, 5000), labels = seq(0, 15, 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Rank in heart gene list (×1,000)", y = NULL) +
  theme_panel +
  theme(legend.position = "none",
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.3),
        plot.margin = margin(t = 0.5, r = 3, b = 0, l = 3, unit = "mm"))

B_heights = c(figX_height_Bcurve_mm, figX_height_Bticks_mm, figX_height_Blollipop_mm)

panel_B = panel_B_curve / panel_B_ticks / panel_B_lollipop +
  plot_layout(heights = B_heights)

ggsave(paste0(fig_folder, "Figure6B.pdf"), panel_B,
       width = figX_width_B_mm / 25.4, height = figX_height_BC_mm / 25.4, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(fig_folder, "Figure6B.png"), panel_B,
       width = figX_width_B_mm / 25.4, height = figX_height_BC_mm / 25.4, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# ---- FIGURE 6C: liver RNA vs protein, jointly significant gene sets by theme ----
#(gene sets with FDR < 0.05 in both liver RNA-seq and liver proteomics, counted
#per theme and coloured by direction RNA / protein)

C_membership = bind_rows(
  tibble(theme = "Lipid, sterol & peroxisome", expected = "Down in both", pathway = c(
    "HALLMARK_FATTY_ACID_METABOLISM", "REACTOME_FATTY_ACID_METABOLISM",
    "REACTOME_PEROXISOMAL_PROTEIN_IMPORT", "REACTOME_ABC_TRANSPORTER_DISORDERS",
    "HALLMARK_PEROXISOME", "HALLMARK_BILE_ACID_METABOLISM",
    "REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT",
    "REACTOME_CHOLESTEROL_BIOSYNTHESIS", "REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS")),
  tibble(theme = "mTORC1", expected = "Down in both", pathway = c(
    "HALLMARK_MTORC1_SIGNALING")),
  tibble(theme = "Redox / hypoxia", expected = "Down in both", pathway = c(
    "REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
    "REACTOME_NUCLEAR_EVENTS_MEDIATED_BY_NFE2L2", "REACTOME_KEAP1_NFE2L2_PATHWAY")),
  tibble(theme = "Immune signalling", expected = "Down in both", pathway = c(
    "REACTOME_CROSS_PRESENTATION_OF_SOLUBLE_EXOGENOUS_ANTIGENS_ENDOSOMES",
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    "REACTOME_DECTIN_1_MEDIATED_NONCANONICAL_NF_KB_SIGNALING",
    "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
    "REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR",
    "REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION", "REACTOME_CLEC7A_DECTIN_1_SIGNALING",
    "REACTOME_C_TYPE_LECTIN_RECEPTORS_CLRS", "REACTOME_INTERLEUKIN_1_SIGNALING")),
  tibble(theme = "Proteasome-linked / cell cycle", expected = "Down in both", pathway = c(
    "REACTOME_PROTEIN_LOCALIZATION",
    "REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT",
    "REACTOME_METABOLISM_OF_POLYAMINES",
    "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",
    "REACTOME_SCF_BETA_TRCP_MEDIATED_DEGRADATION_OF_EMI1",
    "REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA",
    "REACTOME_HEDGEHOG_LIGAND_BIOGENESIS",
    "REACTOME_APC_C_MEDIATED_DEGRADATION_OF_CELL_CYCLE_PROTEINS",
    "REACTOME_DEGRADATION_OF_GLI1_BY_THE_PROTEASOME",
    "REACTOME_DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS", "REACTOME_STABILIZATION_OF_P53",
    "REACTOME_DNA_REPLICATION_PRE_INITIATION",
    "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING",
    "REACTOME_ASYMMETRIC_LOCALIZATION_OF_PCP_PROTEINS", "REACTOME_DEGRADATION_OF_AXIN",
    "REACTOME_REGULATION_OF_RUNX3_EXPRESSION_AND_ACTIVITY",
    "REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS",
    "REACTOME_PCP_CE_PATHWAY", "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS",
    "REACTOME_UCH_PROTEINASES",
    "REACTOME_SWITCHING_OF_ORIGINS_TO_A_POST_REPLICATIVE_STATE",
    "REACTOME_SYNTHESIS_OF_DNA", "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
    "REACTOME_ORC1_REMOVAL_FROM_CHROMATIN",
    "REACTOME_REGULATION_OF_PTEN_STABILITY_AND_ACTIVITY",
    "REACTOME_REGULATION_OF_RUNX2_EXPRESSION_AND_ACTIVITY",
    "REACTOME_REGULATION_OF_RAS_BY_GAPS", "REACTOME_HEDGEHOG_OFF_STATE",
    "REACTOME_SIGNALING_BY_NOTCH4", "REACTOME_SIGNALING_BY_HEDGEHOG",
    "REACTOME_PTEN_REGULATION", "REACTOME_SIGNALING_BY_NOTCH")),
  tibble(theme = "Nuclear pore & SUMOylation", expected = "Up in both", pathway = c(
    "REACTOME_SUMOYLATION",
    "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
    "REACTOME_NUCLEAR_ENVELOPE_BREAKDOWN")),
  tibble(theme = "Cell junctions & actin cytoskeleton", expected = "Up in both", pathway = c(
    "HALLMARK_MITOTIC_SPINDLE", "REACTOME_CELL_CELL_COMMUNICATION",
    "REACTOME_CELL_JUNCTION_ORGANIZATION")),
  tibble(theme = "Extracellular matrix & cell surface", expected = "Up in both", pathway = c(
    "HALLMARK_UV_RESPONSE_DN",
    "REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL")),
  tibble(theme = "Mitochondrial respiration", expected = "Opposite", pathway = c(
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_COMPLEX_I_BIOGENESIS",
    "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM"))
)

C_theme_levels = unique(C_membership$theme)

C_shared = inner_join(
  E_liver %>% select(pathway, NES_rna = NES, padj_rna = padj),
  P_liver %>% select(pathway, NES_prot = NES, padj_prot = padj),
  by = "pathway"
)

C_joint = C_shared %>%
  filter(padj_rna < 0.05, padj_prot < 0.05) %>%
  mutate(direction = case_when(
    NES_rna < 0 & NES_prot < 0 ~ "Down in both",
    NES_rna > 0 & NES_prot > 0 ~ "Up in both",
    TRUE ~ "Opposite"
  ))

stopifnot(nrow(C_shared) == 1104)
stopifnot(all(table(factor(C_joint$direction, levels = names(col_direction))) == c(55, 8, 5)))
stopifnot(!any(duplicated(C_membership$pathway)))
stopifnot(setequal(C_joint$pathway, C_membership$pathway))

C_sets = C_membership %>% left_join(C_joint, by = "pathway")

stopifnot(all(C_sets$direction == C_sets$expected))

C_dat = C_sets %>%
  count(theme, direction, name = "n") %>%
  mutate(theme = factor(theme, levels = C_theme_levels),
         direction = factor(direction, levels = names(col_direction))) %>%
  arrange(theme)

stopifnot(all(C_dat$n == c(9, 1, 4, 9, 32, 3, 3, 2, 5)))

C_rho = cor(C_shared$NES_rna, C_shared$NES_prot, method = "spearman")

stopifnot(round(C_rho, 2) == 0.33)

cat("\n==== Figure 6C: jointly significant liver gene sets by theme ====\n")
print(as.data.frame(C_dat), row.names = FALSE)
cat(sprintf("Spearman rho, liver RNA vs protein NES, all %d shared sets: %.3f\n",
            nrow(C_shared), C_rho))

C_breaks = c(5.5, 8.5)

panel_C = ggplot(C_dat, aes(x = theme, y = n, fill = direction)) +
  geom_vline(xintercept = C_breaks, colour = "grey70", linewidth = 0.3,
             linetype = "dashed") +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2 * 0.75, colour = "grey15") +
  scale_fill_manual(values = col_direction, drop = FALSE,
                    name = "Direction\n(RNA / protein)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = pretty_breaks(n = 4)) +
  labs(x = NULL, y = "Liver significant gene sets\n(RNA + protein)") +
  theme_panel +
  theme(axis.text.x = element_text(size = 6.3, angle = tilt_angle, hjust = 1, vjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification.inside = c(0, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(2.8, "mm"),
        plot.margin = margin(t = 3, r = 3, b = 1, l = 3, unit = "mm"))

ggsave(paste0(fig_folder, "Figure6C.pdf"), panel_C,
       width = figX_width_C_mm / 25.4, height = figX_height_BC_mm / 25.4, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(fig_folder, "Figure6C.png"), panel_C,
       width = figX_width_C_mm / 25.4, height = figX_height_BC_mm / 25.4, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# ---- FIGURE 6A-C: assembled ----
#(panel D, western blots, was added from GraphPad Prism)

fig6_B_stack = panel_B_curve / panel_B_ticks / panel_B_lollipop +
  plot_layout(heights = B_heights)

fig6_bottom = fig6_B_stack |
  wrap_elements(full = panel_C &
                  theme(plot.margin = margin(4, 3, 1, 6, unit = "mm")))

fig6_bottom = fig6_bottom + plot_layout(widths = c(figX_width_B_mm, figX_width_C_mm))

fig6 = (wrap_elements(full = panel_A &
                        theme(plot.margin = margin(t = 5, r = 2, b = 0, l = 2, unit = "mm"))) /
          fig6_bottom) +
  plot_layout(heights = c(figX_height_A_mm, figX_height_BC_mm)) +
  plot_annotation(tag_levels = list(c("A", "B", "", "", "C")),
                  theme = theme(plot.tag = element_text(face = "bold", size = 11),
                                plot.margin = margin(2, 2, 2, 2, unit = "mm")))

ggsave(paste0(fig_folder, "Figure6_ABC.pdf"), fig6,
       width = fig_width_mm / 25.4, height = figX_height_mm / 25.4, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(paste0(fig_folder, "Figure6_ABC.png"), fig6,
       width = fig_width_mm / 25.4, height = figX_height_mm / 25.4, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# ---- SUPPL. FIGURE S4A: liver RNA | protein pathway matrix ----

A7_spec = tribble(
  ~category, ~pathway, ~label, ~ref_rna_nes, ~ref_rna_fdr, ~ref_prot_nes, ~ref_prot_fdr,
  "Lipid & sterol", "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "Cholesterol biosynthesis (R)", -1.88, 8.0e-03, -1.84, 0.040,
  "Lipid & sterol", "HALLMARK_FATTY_ACID_METABOLISM",
  "Fatty-acid metabolism (H)", -2.70, 6.4e-14, -1.73, 0.014,
  "Lipid & sterol", "HALLMARK_PEROXISOME",
  "Peroxisome (H)", -1.93, 2.9e-04, -2.10, 6.9e-04,
  "Lipid & sterol", "HALLMARK_BILE_ACID_METABOLISM",
  "Bile-acid metabolism (H)", -1.88, 9.0e-04, -1.71, 0.033,
  "Lipid & sterol", "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",
  "Activation of gene expression by SREBF/SREBP (R)", -1.10, 0.47, -2.41, 6.9e-05,
  "mTOR", "HALLMARK_MTORC1_SIGNALING",
  "mTORC1 signaling (H)", -1.93, 2.9e-05, -1.98, 1.4e-03,
  "Inflammation", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "TNFα signaling via NF-κB (H)", 2.33, 1.8e-08, 0.98, 0.81
)

A7_dat = bind_rows(
  A7_spec %>% left_join(select(E_liver, pathway, NES, padj, size), by = "pathway") %>%
    mutate(assay = "RNA", ref_nes = ref_rna_nes, ref_fdr = ref_rna_fdr),
  A7_spec %>% left_join(select(P_liver, pathway, NES, padj, size), by = "pathway") %>%
    mutate(assay = "Protein", ref_nes = ref_prot_nes, ref_fdr = ref_prot_fdr)
)

stopifnot(!any(is.na(A7_dat$NES)))
stopifnot(all(abs(round(A7_dat$NES, 2) - A7_dat$ref_nes) < 1e-9 &
                abs(signif(A7_dat$padj, 2) - A7_dat$ref_fdr) < 1e-9 * pmax(A7_dat$ref_fdr, 1e-30)))

A7_dat = A7_dat %>%
  mutate(
    nlp = -log10(padj),
    fdr = factor(ifelse(padj < 0.05, fdr_levels[1], fdr_levels[2]), levels = fdr_levels),
    tissue = "Liver",
    assay = factor(assay, levels = c("RNA", "Protein")),
    category = factor(category, levels = unique(A7_spec$category)),
    label = factor(label, levels = rev(A7_spec$label))
  )

cat("\n==== Suppl. Figure S4A: plotted values ====\n")
print(as.data.frame(A7_dat %>% transmute(category, pathway, assay, NES = round(NES, 4),
                                         FDR = signif(padj, 3), size, fdr)), row.names = FALSE)

panel_S4A = ggplot(A7_dat, aes(x = NES, y = label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_point(data = filter(A7_dat, padj < 0.05),
             aes(size = pmin(nlp, nlp_cap), fill = NES, shape = fdr),
             colour = "grey20", stroke = 0.3) +
  geom_point(data = filter(A7_dat, padj >= 0.05),
             aes(size = pmin(nlp, nlp_cap), colour = NES, shape = fdr),
             fill = "white", stroke = 1.0) +
  scale_fill_gradient2(low = col_nes_low, mid = col_nes_mid, high = col_nes_high,
                       midpoint = 0, limits = nes_limits, oob = squish,
                       breaks = c(-3, 0, 3), name = "NES") +
  scale_colour_gradient2(low = col_nes_low, mid = col_nes_mid, high = col_nes_high,
                         midpoint = 0, limits = nes_limits, oob = squish,
                         guide = "none") +
  scale_size(range = size_range, limits = c(0, nlp_cap), breaks = nlp_breaks,
             labels = c(head(nlp_breaks, -1), paste0("≥", nlp_cap)),
             name = expression(-log[10] * "(FDR)")) +
  scale_shape_manual(values = setNames(c(21, 21), fdr_levels), drop = FALSE,
                     name = NULL) +
  guides(
    fill = guide_colourbar(order = 1, direction = "horizontal",
                           barwidth = unit(16, "mm"),
                           barheight = unit(2.5, "mm"), title.vjust = 0.9),
    size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey70",
                                                       colour = "grey20",
                                                       stroke = 0.3)),
    shape = guide_legend(order = 3, override.aes = list(
      size = 3, fill = c("grey55", "white"), colour = c("grey20", "grey40"),
      stroke = c(0.3, 1.0)))
  ) +
  scale_x_continuous(limits = nes_limits, breaks = seq(-3, 3, 1.5),
                     labels = c("-3", "", "0", "", "3")) +
  scale_y_discrete(position = "right") +
  facet_nested(tissue + category ~ assay, scales = "free_y", space = "free_y",
               switch = "y",
               nest_line = element_line(colour = "grey30", linewidth = 0.4),
               strip = strip_nested(
                 text_y = elem_list_text(face = c("bold", "plain"),
                                         size = c(8, 7),
                                         angle = c(0, 0), hjust = c(1, 1)),
                 by_layer_y = TRUE)) +
  coord_cartesian(clip = "off") +
  labs(x = "NES", y = NULL) +
  theme_panel +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.x = unit(1.5, "mm"),
        legend.spacing.y = unit(0.5, "mm"),
        legend.margin = margin(0, 2, 0, 0, unit = "mm"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
        axis.text.y.right = element_text(size = 6.8, hjust = 0),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  theme(strip.text.x = element_text(size = 8, face = "bold"),
        panel.spacing.x = unit(2.5, "mm"),
        panel.spacing.y = unit(2.2, "mm"),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.3),
        legend.position = "right",
        legend.justification = c(0, 0.5),
        legend.box.margin = margin(0, 0, 0, 4, unit = "mm"))

ggsave(paste0(fig_folder, "SupplFigureS4A.png"), panel_S4A,
       width = fig_width_mm / 25.4, height = fig7_height_A_mm / 25.4, units = "in",
       dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                    ---- SUPPLEMENTARY TABLES S1-S3 ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

hs = createStyle(textDecoration = "bold", fgFill = "#E7EEF0", border = "bottom")

enrich_tables = bind_rows(
  E_cortex %>% mutate(Tissue = "Cortex"),
  E_heart %>% mutate(Tissue = "Heart"),
  E_liver %>% mutate(Tissue = "Liver"),
  P_liver %>% mutate(Tissue = "Liver protein")
) %>%
  transmute(Tissue,
            Collection = ifelse(grepl("^HALLMARK_", pathway), "Hallmark",
                                ifelse(grepl("^REACTOME_", pathway), "Reactome", NA)),
            `Gene set` = sub("^(HALLMARK|REACTOME)_", "", pathway),
            pathway,
            NES,
            ES,
            `P value` = pval,
            FDR = padj,
            `Set size` = size,
            `Leading-edge genes` = leadingEdge)

# ---- Table S1: RNA-seq gene-set enrichment ----

#(table legends are in the supplementary file of the manuscript)

wb = createWorkbook()

for (tissue in c("Cortex", "Heart", "Liver")) {

  S1 = enrich_tables %>%
    filter(Tissue == tissue) %>%
    select(-Tissue, -pathway) %>%
    arrange(FDR, desc(abs(NES)))

  addWorksheet(wb, tissue)
  writeData(wb, tissue, S1, headerStyle = hs)
  freezePane(wb, tissue, firstRow = TRUE)
  setColWidths(wb, tissue, cols = seq_along(S1), widths = c(10, 60, 9, 9, 11, 11, 8, 100))

}

saveWorkbook(wb, paste0(tab_folder, "Suppl Table S1 - RNA-seq gene-set enrichment.xlsx"), overwrite = TRUE)

# ---- Table S2: complex I biogenesis, heart RNA and protein ----

ci = sort(set_genes)

hr = DEGS[DEGS$Region == "Hr", ]
hr = hr[!duplicated(hr$gene_symbol), ]

prh = DAPS_HT[DAPS_HT$ident.1 == "4W", ]
prh = prh[!duplicated(prh$GeneSymbol), ]

r_i = match(ci, hr$gene_symbol)
p_i = match(ci, prh$GeneSymbol)

S2 = data.frame(Gene = ci,
                `NDUF subunit or assembly factor (Nduf*)` = ifelse(grepl("^Nduf", ci), "yes", "no"),
                `RNA log2FC` = hr$logFC[r_i],
                `RNA P value` = hr$P.Value[r_i],
                `RNA FDR` = hr$adj.P.Val[r_i],
                `RNA leading edge (Figure 6B)` = ifelse(ci %in% B_le, "yes", "no"),
                `Protein UniProt ID` = prh$name[p_i],
                `Protein log2FC` = prh$diff[p_i],
                `Protein P value` = prh$pval[p_i],
                `Protein FDR` = prh$adj_pval[p_i],
                check.names = FALSE)

S2 = S2[order(S2$`RNA leading edge (Figure 6B)` != "yes", -S2$`RNA log2FC`), ]

S2_le = S2[S2$`RNA leading edge (Figure 6B)` == "yes", ]
S2_prot = S2[!is.na(S2$`Protein log2FC`), ]
S2_nduf = S2_prot[S2_prot$`NDUF subunit or assembly factor (Nduf*)` == "yes", ]

cat(sprintf("\nS2 check: leading edge %d (all up: %s; FDR<0.05: %d); proteins detected %d (Nduf* %d, up %d, median %.3f; min protein P %.3f; protein FDR<0.05: %d)\n",
            nrow(S2_le), all(S2_le$`RNA log2FC` > 0), sum(S2_le$`RNA FDR` < 0.05), nrow(S2_prot),
            nrow(S2_nduf), sum(S2_nduf$`Protein log2FC` > 0), median(S2_nduf$`Protein log2FC`),
            min(S2_prot$`Protein P value`), sum(S2_prot$`Protein FDR` < 0.05)))

wb = createWorkbook()

addWorksheet(wb, "Complex I biogenesis")
writeData(wb, "Complex I biogenesis", S2, headerStyle = hs)
freezePane(wb, "Complex I biogenesis", firstRow = TRUE)
setColWidths(wb, "Complex I biogenesis", cols = seq_along(S2), widths = c(12, 18, 11, 12, 12, 14, 14, 13, 13, 13))

saveWorkbook(wb, paste0(tab_folder, "Suppl Table S2 - Heart complex I RNA and protein.xlsx"), overwrite = TRUE)

# ---- Table S3: liver proteome enrichment and RNA vs protein ----

S3_prot = enrich_tables %>% filter(Tissue == "Liver protein") %>% select(-Tissue)
S3_rna  = enrich_tables %>% filter(Tissue == "Liver") %>% select(-Tissue)

S3_shared = intersect(S3_rna$pathway, S3_prot$pathway)
S3_a = S3_rna[match(S3_shared, S3_rna$pathway), ]
S3_b = S3_prot[match(S3_shared, S3_prot$pathway), ]

S3_cmp = data.frame(Collection = S3_a$Collection,
                    `Gene set` = S3_a$`Gene set`,
                    `RNA NES` = S3_a$NES,
                    `RNA FDR` = S3_a$FDR,
                    `Protein NES` = S3_b$NES,
                    `Protein FDR` = S3_b$FDR,
                    check.names = FALSE)

S3_both = S3_cmp$`RNA FDR` < 0.05 & S3_cmp$`Protein FDR` < 0.05
S3_cmp$`FDR < 0.05 in both` = ifelse(S3_both, "yes", "no")
S3_cmp$Direction = ifelse(!S3_both, NA_character_, ifelse(sign(S3_cmp$`RNA NES`) == sign(S3_cmp$`Protein NES`), "same", "opposite"))

cat(sprintf("S3 check: shared %d, both significant %d, same direction %d, opposite %d, Spearman rho %.2f\n",
            nrow(S3_cmp), sum(S3_both), sum(S3_cmp$Direction == "same", na.rm = TRUE),
            sum(S3_cmp$Direction == "opposite", na.rm = TRUE),
            cor(S3_cmp$`RNA NES`, S3_cmp$`Protein NES`, method = "spearman")))

S3_cmp = S3_cmp[order(S3_cmp$`FDR < 0.05 in both` != "yes", S3_cmp$`RNA FDR`), ]

srebp_le = strsplit(S3_prot$`Leading-edge genes`[S3_prot$pathway == "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP"], ", ")[[1]]

prl = DAPS_LV[DAPS_LV$ident.1 == "4W", ]
prl = prl[!duplicated(prl$GeneSymbol), ]

lv = DEGS[DEGS$Region == "Lv", ]
lv = lv[!duplicated(lv$gene_symbol), ]

p_i = match(srebp_le, prl$GeneSymbol)
r_i = match(srebp_le, lv$gene_symbol)

S3_sre = data.frame(Gene = srebp_le,
                    `Protein UniProt ID` = prl$name[p_i],
                    `Protein log2FC` = prl$diff[p_i],
                    `Protein P value` = prl$pval[p_i],
                    `Protein FDR` = prl$adj_pval[p_i],
                    `RNA log2FC` = lv$logFC[r_i],
                    `RNA P value` = lv$P.Value[r_i],
                    `RNA FDR` = lv$adj.P.Val[r_i],
                    check.names = FALSE)

cat("SREBP leading edge:", nrow(S3_sre), "proteins; protein FDR<0.05:", sum(S3_sre$`Protein FDR` < 0.05),
    "; all protein log2FC < 0:", all(S3_sre$`Protein log2FC` < 0), "\n")

wb = createWorkbook()

S3_prot_sheet = S3_prot %>% select(-pathway) %>% arrange(FDR, desc(abs(NES)))

addWorksheet(wb, "Protein enrichment")
writeData(wb, "Protein enrichment", S3_prot_sheet, headerStyle = hs)
freezePane(wb, "Protein enrichment", firstRow = TRUE)
setColWidths(wb, "Protein enrichment", cols = seq_along(S3_prot_sheet), widths = c(10, 60, 9, 9, 11, 11, 8, 100))

addWorksheet(wb, "RNA vs protein")
writeData(wb, "RNA vs protein", S3_cmp, headerStyle = hs)
freezePane(wb, "RNA vs protein", firstRow = TRUE)
setColWidths(wb, "RNA vs protein", cols = seq_along(S3_cmp), widths = c(10, 60, 9, 11, 11, 11, 16, 10))

addWorksheet(wb, "SREBP-target proteins")
writeData(wb, "SREBP-target proteins", S3_sre, headerStyle = hs)
freezePane(wb, "SREBP-target proteins", firstRow = TRUE)
setColWidths(wb, "SREBP-target proteins", cols = seq_along(S3_sre), widths = c(10, 16, 13, 13, 12, 11, 12, 12))

saveWorkbook(wb, paste0(tab_folder, "Suppl Table S3 - Liver proteome and RNA comparison.xlsx"), overwrite = TRUE)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#                    ---- SUPPLEMENTARY FIGURE S4 ----
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# ---- Panel B: SREBF/SREBP leading-edge proteins, protein and RNA log2FC ----

S4_long = rbind(
  data.frame(Gene = S3_sre$Gene, Layer = "Protein", lfc = S3_sre$`Protein log2FC`, fdr = S3_sre$`Protein FDR`),
  data.frame(Gene = S3_sre$Gene, Layer = "RNA", lfc = S3_sre$`RNA log2FC`, fdr = S3_sre$`RNA FDR`))
S4_long = S4_long[!is.na(S4_long$lfc), ]
S4_long$Gene = factor(S4_long$Gene, levels = rev(S3_sre$Gene[order(S3_sre$`Protein log2FC`)]))
S4_long$sig = ifelse(S4_long$fdr < 0.05, "FDR < 0.05", "FDR ≥ 0.05")

panel_S4B = ggplot(S4_long, aes(lfc, Gene, colour = Layer)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(fill = interaction(Layer, sig), shape = Layer), size = 3, stroke = 1,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Protein" = "#1b7f86", "RNA" = "#8a8a8a"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_fill_manual(values = c("Protein.FDR < 0.05" = "#1b7f86", "RNA.FDR < 0.05" = "#8a8a8a",
                               "Protein.FDR ≥ 0.05" = "white", "RNA.FDR ≥ 0.05" = "white"),
                    guide = "none") +
  labs(x = expression(log[2]~"fold change (noise vs control)"), y = NULL,
       subtitle = "Leading-edge proteins: Activation of gene expression by SREBF/SREBP (R)") +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "right",
        plot.subtitle = element_text(size = 10))

# ---- Assembled (panel A is the 600 dpi export of the liver RNA | protein matrix) ----

panel_S4A_img = wrap_elements(full = rasterGrob(readPNG(paste0(fig_folder, "SupplFigureS4A.png")),
                                                interpolate = TRUE))

fig_S4 = panel_S4A_img / panel_S4B +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = "Arial"))

agg_png(paste0(fig_folder, "Suppl_Figure_S4.png"), width = 7.2, height = 6.6, units = "in", res = 300)
print(fig_S4)
invisible(dev.off())

cat("\nDone. Results written to: ", normalizePath(res_folder), "\n", sep = "")
