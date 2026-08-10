# ============================================================================
# DynaMac-GP v2.3 — direct-input, nominal-P tiering version
#
# Directly reads:
#   1) one1k_finngen_R12_TEMPOROMANDIB_INCLAVO.csv              [all MR]
#   2) macro_finngen_R12_TEMPOROMANDIB_INCLAVO.csv              [all MR]
#   3) one1k_coloc_summary_coloc_abf.csv                         [coloc]
#   4) MacroMap_batch_coloc_summary.tsv                          [coloc]
#
# Optional QC-only files:
#   5) one1k_finngen_R12_TEMPOROMANDIB_INCLAVO_pval_mr_less_0.05.csv
#   6) macro_finngen_R12_TEMPOROMANDIB_INCLAVO_pval_mr_less_0.05.csv
#

# ============================================================================

# =========================
# 0. Packages
# =========================
cran_pkgs <- c(
  "tidyverse", "readr", "openxlsx", "pheatmap", "cluster",
  "ggrepel", "scales", "patchwork"
)

required_bioc <- c(
  "ConsensusClusterPlus", "AnnotationDbi", "org.Hs.eg.db"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in required_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(openxlsx)
  library(pheatmap)
  library(cluster)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(ConsensusClusterPlus)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# Explicitly prefer dplyr verbs after loading Bioconductor packages.
# This prevents AnnotationDbi::select() and IRanges::slice() from being
# dispatched on tibble objects.
select <- dplyr::select
slice <- dplyr::slice
filter <- dplyr::filter
rename <- dplyr::rename
rename_with <- dplyr::rename_with
mutate <- dplyr::mutate
transmute <- dplyr::transmute
summarise <- dplyr::summarise
arrange <- dplyr::arrange
group_by <- dplyr::group_by
group_modify <- dplyr::group_modify
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
right_join <- dplyr::right_join
full_join <- dplyr::full_join
inner_join <- dplyr::inner_join
anti_join <- dplyr::anti_join
distinct <- dplyr::distinct
count <- dplyr::count
pull <- dplyr::pull
case_when <- dplyr::case_when
if_else <- dplyr::if_else

# Optional enrichment packages. Set FALSE if package installation is undesired.
run_enrichment <- TRUE
if (run_enrichment) {
  optional_bioc <- c("clusterProfiler", "enrichplot")
  for (p in optional_bioc) {
    if (!requireNamespace(p, quietly = TRUE)) {
      try(BiocManager::install(p, ask = FALSE, update = FALSE), silent = TRUE)
    }
  }
}

# =========================
# 1. Paths — only edit input_dir and out_dir
# =========================
input_dir <- "G:/hz/1.11汇总/ly/tmdv1/"
out_dir   <- "G:/hz/1.11汇总/ly/tmdv1/re2/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  tables      = file.path(out_dir, "01_tables"),
  qc          = file.path(out_dir, "02_QC_figures"),
  main        = file.path(out_dir, "03_main_clustering"),
  tier        = file.path(out_dir, "04_tier_results"),
  sensitivity = file.path(out_dir, "05_sensitivity"),
  enrichment  = file.path(out_dir, "06_enrichment")
)
walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

files <- list(
  onek1k_mr = file.path(
    input_dir,
    "o_table_s2_mr_results.csv"
  ),
  macromap_mr = file.path(
    input_dir,
    "m_table_s2_mr_results.csv"
  ),
  onek1k_coloc = file.path(
    input_dir,
    "onek1k_coloc_summary_coloc_abf.csv"
  ),
  macromap_coloc = file.path(
    input_dir,
    "MacroMap_batch_coloc_summary.tsv"
  ),
  onek1k_selected_qc = file.path(
    input_dir,
    "o_0.05table_s2_mr_results.csv"
  ),
  macromap_selected_qc = file.path(
    input_dir,
    "m_0.05table_s2_mr_results.csv"
  )
)

# =========================
# 2. Main parameters
# =========================
set.seed(12345)

mr_screen_p <- 0.05
min_p <- 1e-300
max_logp <- 50
strong_coloc_cutoff <- 0.70
moderate_coloc_cutoff <- 0.50

# How to choose one representative MR row when a gene-context has multiple
# methods in OneK1K:
#   "coloc_compatible": if any method has P<0.05, select the highest-priority
#                       method among those passing P<0.05. This matches the
#                       actual uploaded coloc set (713 unique pairs).
#   "primary_method": always prioritize Wald ratio for 1 SNP and IVW for >1
#                     SNP. Six uploaded OneK1K coloc pairs then become
#                     non-primary-method screens and are reported as QC.
mr_method_selection_mode <- "coloc_compatible"

candidate_mode <- "common_any_coloc"
main_aggregation <- "top_quartile"
top_fraction <- 0.25
min_context_coverage_within_module <- 0.25
min_module_coverage_per_block <- 0.50
main_static_weight <- 0.50
main_dynamic_weight <- 0.50
main_distance <- "euclidean"
main_reps <- 1000
sensitivity_reps <- 500
p_item <- 0.80
p_feature <- 0.80
max_k_cap <- 6
min_cluster_size <- 2
tier_mr_p <- 0.05

# =========================
# 3. Utility functions
# =========================
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_max <- function(x) {
  x <- safe_numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- safe_numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

rank_percentile <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n_ok <- sum(ok)
  if (n_ok == 0) return(out)
  out[ok] <- (rank(x[ok], ties.method = "average") - 0.5) / n_ok
  out
}

top_fraction_mean <- function(x, fraction = 0.25) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  k <- max(1L, ceiling(length(x) * fraction))
  mean(sort(x, decreasing = TRUE)[seq_len(k)])
}

top_fraction_signed <- function(x, fraction = 0.25) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  k <- max(1L, ceiling(length(x) * fraction))
  idx <- order(abs(x), decreasing = TRUE)[seq_len(k)]
  mean(x[idx])
}

aggregate_vector <- function(x, method = "top_quartile", signed = FALSE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (method == "top_quartile") {
    if (signed) return(top_fraction_signed(x, top_fraction))
    return(top_fraction_mean(x, top_fraction))
  }
  if (method == "max") {
    if (signed) return(x[which.max(abs(x))])
    return(max(x))
  }
  if (method == "median") return(median(x))
  stop("Unknown aggregation method: ", method)
}

adjusted_rand_index <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.factor(x[ok]); y <- as.factor(y[ok])
  n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y)
  comb2 <- function(z) z * (z - 1) / 2
  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_ai * sum_bj / total
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (denom == 0) return(NA_real_)
  (sum_nij - expected) / denom
}

normalized_mutual_information <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.factor(x[ok]); y <- as.factor(y[ok])
  tab <- table(x, y)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(pxy[nz] * log(pxy[nz] / (px[nz[, 1]] * py[nz[, 2]])))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

coclustering_jaccard <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2) return(NA_real_)
  pairs <- combn(seq_len(n), 2)
  a <- x[pairs[1, ]] == x[pairs[2, ]]
  b <- y[pairs[1, ]] == y[pairs[2, ]]
  union_n <- sum(a | b)
  if (union_n == 0) return(NA_real_)
  sum(a & b) / union_n
}

make_dist <- function(mat, method) {
  if (method == "euclidean") return(dist(mat, method = "euclidean"))
  if (method == "pearson") {
    r <- cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
    r[is.na(r)] <- 0
    return(as.dist(1 - r))
  }
  stop("Unsupported distance: ", method)
}

calc_pac <- function(consensus_matrix, lower = 0.1, upper = 0.9) {
  cm <- consensus_matrix[upper.tri(consensus_matrix)]
  mean(cm > lower & cm < upper)
}

save_plot_both <- function(plot_obj, filename_no_ext, width = 8, height = 6) {
  ggsave(paste0(filename_no_ext, ".png"), plot_obj,
         width = width, height = height, dpi = 300)
  ggsave(paste0(filename_no_ext, ".pdf"), plot_obj,
         width = width, height = height)
}

mr_method_rank <- function(method, nsnp) {
  method_lower <- tolower(as.character(method))
  case_when(
    !is.na(nsnp) & nsnp <= 1 & str_detect(method_lower, "wald") ~ 1L,
    !is.na(nsnp) & nsnp > 1 & str_detect(method_lower, "inverse|ivw") ~ 1L,
    str_detect(method_lower, "wald") ~ 2L,
    str_detect(method_lower, "inverse|ivw") ~ 2L,
    str_detect(method_lower, "maximum|likelihood") ~ 3L,
    str_detect(method_lower, "weighted") ~ 4L,
    str_detect(method_lower, "egger") ~ 5L,
    TRUE ~ 9L
  )
}

select_representative_mr <- function(df, mode = "coloc_compatible") {
  df2 <- df %>%
    mutate(
      method_rank = mr_method_rank(method, nsnp),
      any_method_selected = any(pval < mr_screen_p, na.rm = TRUE),
      primary_candidate = method_rank == min(method_rank, na.rm = TRUE)
    )

  if (mode == "coloc_compatible" && any(df2$pval < mr_screen_p, na.rm = TRUE)) {
    df2 <- df2 %>% filter(pval < mr_screen_p)
  }

  df2 %>%
    arrange(method_rank, pval) %>%
    slice(1) %>%
    select(-method_rank, -primary_candidate)
}

select_all_representative_mr <- function(df, mode) {
  df %>%
    group_by(source, gene_symbol, context) %>%
    group_modify(~select_representative_mr(.x, mode = mode)) %>%
    ungroup()
}

map_ensembl_to_symbol <- function(ensembl_base) {
  keys <- unique(ensembl_base[!is.na(ensembl_base) & ensembl_base != ""])
  if (length(keys) == 0) return(character())
  suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = keys,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
  )
}

canonicalize_gene_symbols <- function(symbols) {
  keys <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(keys) == 0) return(symbols)

  direct <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = keys,
      keytype = "SYMBOL",
      column = "SYMBOL",
      multiVals = "first"
    )
  )

  missing_keys <- keys[is.na(direct[keys])]
  alias <- setNames(rep(NA_character_, length(missing_keys)), missing_keys)
  if (length(missing_keys) > 0) {
    alias <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = missing_keys,
        keytype = "ALIAS",
        column = "SYMBOL",
        multiVals = "first"
      )
    )
  }

  lookup <- setNames(keys, keys)
  lookup[names(direct)[!is.na(direct)]] <- unname(direct[!is.na(direct)])
  lookup[names(alias)[!is.na(alias)]] <- unname(alias[!is.na(alias)])
  unname(lookup[symbols])
}

# =========================
# 4. Read, parse and harmonize the uploaded raw files
# =========================
required_files <- files[c("onek1k_mr", "macromap_mr", "onek1k_coloc", "macromap_coloc")]
missing_required <- names(required_files)[!file.exists(unlist(required_files))]
if (length(missing_required) > 0) {
  stop(
    "Missing required file(s): ", paste(missing_required, collapse = ", "),
    "\nExpected paths:\n", paste(unlist(required_files), collapse = "\n")
  )
}

onek1k_mr_raw <- read_csv(files$onek1k_mr, show_col_types = FALSE)
macromap_mr_raw <- read_csv(files$macromap_mr, show_col_types = FALSE)
onek1k_coloc_raw <- read_csv(files$onek1k_coloc, show_col_types = FALSE)
macromap_coloc_raw <- read_tsv(files$macromap_coloc, show_col_types = FALSE)

# ---- OneK1K all-MR table ----
onek1k_mr_long <- onek1k_mr_raw %>%
  transmute(
    source = "OneK1K_static",
    exposure_id = as.character(id.exposure),
    context = str_extract(exposure_id, "^[^_]+"),
    gene_symbol_original = str_replace(exposure_id, "^[^_]+_", ""),
    gene_symbol = canonicalize_gene_symbols(gene_symbol_original),
    gene_id = NA_character_,
    method = as.character(method),
    nsnp = safe_numeric(nsnp),
    beta = safe_numeric(b),
    se = safe_numeric(se_mr),
    pval = safe_numeric(pval_mr)
  ) %>%
  filter(
    !is.na(gene_symbol), gene_symbol != "",
    !is.na(context), context != "",
    is.finite(beta), is.finite(se), is.finite(pval),
    pval >= 0, pval <= 1
  )

# ---- MacroMap all-MR table ----
macromap_mr_parsed <- macromap_mr_raw %>%
  transmute(
    source = "Macrophage_dynamic",
    exposure_id = as.character(id.exposure),
    context = str_extract(exposure_id, "^[^:]+") %>%
      str_remove("\\.permuted$"),
    gene_id_versioned = str_replace(exposure_id, "^.*::", ""),
    gene_id = str_remove(gene_id_versioned, "\\..*$"),
    method = as.character(method),
    nsnp = safe_numeric(nsnp),
    beta = safe_numeric(b),
    se = safe_numeric(se),
    pval = safe_numeric(pval_mr)
  ) %>%
  filter(
    !is.na(gene_id), gene_id != "",
    !is.na(context), context != "",
    is.finite(beta), is.finite(se), is.finite(pval),
    pval >= 0, pval <= 1
  )

ensembl_map <- map_ensembl_to_symbol(macromap_mr_parsed$gene_id)

macro_gene_map <- macromap_mr_parsed %>%
  distinct(gene_id, gene_id_versioned) %>%
  mutate(
    mapped_symbol = unname(ensembl_map[gene_id]),
    gene_symbol = if_else(
      !is.na(mapped_symbol) & mapped_symbol != "",
      mapped_symbol,
      gene_id
    ),
    mapping_status = if_else(
      gene_symbol == gene_id,
      "Unmapped_Ensembl_retained",
      "Mapped_to_symbol"
    )
  )

write_csv(macro_gene_map, file.path(dirs$tables, "MacroMap_Ensembl_to_symbol_mapping.csv"))

macromap_mr_long <- macromap_mr_parsed %>%
  left_join(
    macro_gene_map %>% select(gene_id, gene_id_versioned, gene_symbol, mapping_status),
    by = c("gene_id", "gene_id_versioned")
  ) %>%
  select(
    source, gene_symbol, context, beta, se, pval, method, nsnp,
    exposure_id, gene_id, gene_id_versioned, mapping_status
  )

# Compare fixed-primary and coloc-compatible method selection in OneK1K.
onek1k_primary_rows <- onek1k_mr_long %>%
  select_all_representative_mr(mode = "primary_method") %>%
  select(
    gene_symbol, context,
    primary_method = method, primary_nsnp = nsnp,
    primary_beta = beta, primary_se = se, primary_pval = pval
  )

onek1k_compatible_rows <- onek1k_mr_long %>%
  select_all_representative_mr(mode = "coloc_compatible") %>%
  select(
    gene_symbol, context, any_method_selected,
    compatible_method = method, compatible_nsnp = nsnp,
    compatible_beta = beta, compatible_se = se, compatible_pval = pval
  )

method_selection_qc <- full_join(
  onek1k_primary_rows, onek1k_compatible_rows,
  by = c("gene_symbol", "context")
) %>%
  mutate(
    primary_selected = primary_pval < mr_screen_p,
    compatible_selected = compatible_pval < mr_screen_p,
    selection_discordant = primary_selected != compatible_selected |
      primary_method != compatible_method
  )

write_csv(
  method_selection_qc %>% filter(selection_discordant),
  file.path(dirs$tables, "OneK1K_MR_method_selection_discordance.csv")
)

# Choose one representative MR method per resource-gene-context.
mr_all <- bind_rows(onek1k_mr_long, macromap_mr_long) %>%
  select_all_representative_mr(mode = mr_method_selection_mode)

# ---- OneK1K coloc table ----
one_match <- str_match(
  as.character(onek1k_coloc_raw$gwas1),
  "^onek1k_([^_]+)_(.+?)_coloc_"
)

onek1k_coloc <- onek1k_coloc_raw %>%
  mutate(
    source = "OneK1K_static",
    context = one_match[, 2],
    gene_symbol_original = one_match[, 3],
    gene_symbol = canonicalize_gene_symbols(gene_symbol_original),
    exposure_id = as.character(gwas1),
    PP.H3 = safe_numeric(PP.H3.abf),
    PP.H4 = safe_numeric(PP.H4.abf),
    lead_snp = NA_character_,
    n_coloc_snps = safe_numeric(nsnps)
  ) %>%
  filter(
    !is.na(gene_symbol), gene_symbol != "",
    !is.na(context), context != "",
    is.finite(PP.H4)
  ) %>%
  select(
    source, gene_symbol, context, PP.H3, PP.H4,
    exposure_id, lead_snp, n_coloc_snps
  )

# ---- MacroMap coloc table ----
macromap_coloc <- macromap_coloc_raw %>%
  transmute(
    source = "Macrophage_dynamic",
    context = str_remove(as.character(condition_short), "\\.permuted$"),
    gene_id_versioned = as.character(gene_id),
    gene_id = str_remove(gene_id_versioned, "\\..*$"),
    PP.H3 = safe_numeric(PP.H3.abf),
    PP.H4 = safe_numeric(PP.H4.abf),
    exposure_id = paste(context, gene_id_versioned, sep = "::"),
    lead_snp = as.character(top_eqtl_snp),
    n_coloc_snps = safe_numeric(nsnps),
    status = as.character(status)
  ) %>%
  filter(is.na(status) | status == "ok") %>%
  left_join(
    macro_gene_map %>% distinct(gene_id, gene_symbol),
    by = "gene_id"
  ) %>%
  mutate(
    gene_symbol = if_else(
      is.na(gene_symbol) | gene_symbol == "",
      gene_id,
      gene_symbol
    )
  ) %>%
  filter(
    !is.na(gene_symbol), gene_symbol != "",
    !is.na(context), context != "",
    is.finite(PP.H4)
  ) %>%
  select(
    source, gene_symbol, context, PP.H3, PP.H4,
    exposure_id, lead_snp, n_coloc_snps
  )

coloc_all <- bind_rows(onek1k_coloc, macromap_coloc) %>%
  group_by(source, gene_symbol, context) %>%
  arrange(desc(PP.H4), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# ---- Context-to-module mapping generated from the actual uploaded contexts ----
context_map <- tribble(
  ~source, ~context, ~module, ~block,
  "OneK1K_static", "bin",       "Static_B_cell",          "Static",
  "OneK1K_static", "bmem",      "Static_B_cell",          "Static",
  "OneK1K_static", "plasma",    "Static_B_cell",          "Static",
  "OneK1K_static", "cd4nc",     "Static_CD4_T",           "Static",
  "OneK1K_static", "cd4et",     "Static_CD4_T",           "Static",
  "OneK1K_static", "cd4sox4",   "Static_CD4_T",           "Static",
  "OneK1K_static", "cd8nc",     "Static_CD8_T",           "Static",
  "OneK1K_static", "cd8et",     "Static_CD8_T",           "Static",
  "OneK1K_static", "cd8s100b",  "Static_CD8_T",           "Static",
  "OneK1K_static", "nk",        "Static_NK",              "Static",
  "OneK1K_static", "nkr",       "Static_NK",              "Static",
  "OneK1K_static", "monoc",     "Static_Monocyte",        "Static",
  "OneK1K_static", "mononc",    "Static_Monocyte",        "Static",
  "OneK1K_static", "dc",        "Static_Dendritic_cell",  "Static",
  "Macrophage_dynamic", "Ctrl_6",    "Macrophage_Ctrl",                "Dynamic",
  "Macrophage_dynamic", "Ctrl_24",   "Macrophage_Ctrl",                "Dynamic",
  "Macrophage_dynamic", "IFNG_6",    "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "IFNG_24",   "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "sLPS_6",    "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "sLPS_24",   "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "P3C_6",     "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "P3C_24",    "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "CIL_6",     "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "CIL_24",    "Macrophage_ProInflam",           "Dynamic",
  "Macrophage_dynamic", "IL4_6",     "Macrophage_Repair_Resolution",   "Dynamic",
  "Macrophage_dynamic", "IL4_24",    "Macrophage_Repair_Resolution",   "Dynamic",
  "Macrophage_dynamic", "LIL10_6",   "Macrophage_Repair_Resolution",   "Dynamic",
  "Macrophage_dynamic", "LIL10_24",  "Macrophage_Repair_Resolution",   "Dynamic",
  "Macrophage_dynamic", "IFNB_6",    "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "IFNB_24",   "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "R848_6",    "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "R848_24",   "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "PIC_6",     "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "PIC_24",    "Macrophage_IFN_Antiviral",      "Dynamic",
  "Macrophage_dynamic", "MBP_6",     "Macrophage_NeuroLike",          "Dynamic",
  "Macrophage_dynamic", "MBP_24",    "Macrophage_NeuroLike",          "Dynamic",
  "Macrophage_dynamic", "Prec_D0",   "Macrophage_Precursor",          "Dynamic",
  "Macrophage_dynamic", "Prec_D2",   "Macrophage_Precursor",          "Dynamic"
) %>%
  mutate(
    module = make.names(module, unique = FALSE),
    feature = paste(block, module, sep = "__")
  )

write_csv(context_map, file.path(dirs$tables, "context_module_map_used.csv"))

# Validate all observed contexts.
unmapped_contexts <- mr_all %>%
  distinct(source, context) %>%
  anti_join(context_map, by = c("source", "context"))
if (nrow(unmapped_contexts) > 0) {
  write_csv(unmapped_contexts, file.path(dirs$tables, "ERROR_unmapped_contexts.csv"))
  stop("Some MR contexts are absent from the internally defined context map.")
}

# Coloc rows that do not match a representative MR row are retained for QC.
coloc_without_mr <- coloc_all %>%
  anti_join(mr_all, by = c("source", "gene_symbol", "context"))
write_csv(
  coloc_without_mr,
  file.path(dirs$tables, "coloc_rows_without_matching_representative_MR.csv")
)

# Optional uploaded P<0.05 files are used only to confirm preprocessing.
qc_selected_summary <- tibble()
if (file.exists(files$onek1k_selected_qc)) {
  one_sel_raw <- read_csv(files$onek1k_selected_qc, show_col_types = FALSE)
  one_sel <- one_sel_raw %>%
    transmute(
      source = "OneK1K_static",
      exposure_id = as.character(id.exposure),
      context = str_extract(exposure_id, "^[^_]+"),
      gene_symbol = canonicalize_gene_symbols(
        str_replace(exposure_id, "^[^_]+_", "")
      )
    ) %>%
    distinct(source, gene_symbol, context)
  qc_selected_summary <- bind_rows(
    qc_selected_summary,
    tibble(
      resource = "OneK1K",
      selected_file_unique_pairs = nrow(one_sel),
      coloc_unique_pairs = nrow(onek1k_coloc),
      selected_without_coloc = nrow(anti_join(
        one_sel, onek1k_coloc,
        by = c("source", "gene_symbol", "context")
      )),
      coloc_without_selected = nrow(anti_join(
        onek1k_coloc, one_sel,
        by = c("source", "gene_symbol", "context")
      ))
    )
  )
}

if (file.exists(files$macromap_selected_qc)) {
  macro_sel_raw <- read_csv(files$macromap_selected_qc, show_col_types = FALSE)
  macro_sel <- macro_sel_raw %>%
    transmute(
      source = "Macrophage_dynamic",
      exposure_id = as.character(id.exposure),
      context = str_extract(exposure_id, "^[^:]+") %>%
        str_remove("\\.permuted$"),
      gene_id = str_replace(exposure_id, "^.*::", "") %>%
        str_remove("\\..*$")
    ) %>%
    left_join(
      macro_gene_map %>% distinct(gene_id, gene_symbol),
      by = "gene_id"
    ) %>%
    mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol)) %>%
    distinct(source, gene_symbol, context)
  qc_selected_summary <- bind_rows(
    qc_selected_summary,
    tibble(
      resource = "MacroMap",
      selected_file_unique_pairs = nrow(macro_sel),
      coloc_unique_pairs = nrow(macromap_coloc),
      selected_without_coloc = nrow(anti_join(
        macro_sel, macromap_coloc,
        by = c("source", "gene_symbol", "context")
      )),
      coloc_without_selected = nrow(anti_join(
        macromap_coloc, macro_sel,
        by = c("source", "gene_symbol", "context")
      ))
    )
  )
}

if (nrow(qc_selected_summary) > 0) {
  write_csv(
    qc_selected_summary,
    file.path(dirs$tables, "uploaded_selected_files_QC.csv")
  )
}

# Input audit summary.
input_audit <- bind_rows(
  tibble(
    resource = "OneK1K",
    raw_MR_rows = nrow(onek1k_mr_raw),
    representative_MR_pairs = sum(mr_all$source == "OneK1K_static"),
    coloc_pairs = nrow(onek1k_coloc),
    contexts = n_distinct(onek1k_mr_long$context),
    unique_genes = n_distinct(onek1k_mr_long$gene_symbol)
  ),
  tibble(
    resource = "MacroMap",
    raw_MR_rows = nrow(macromap_mr_raw),
    representative_MR_pairs = sum(mr_all$source == "Macrophage_dynamic"),
    coloc_pairs = nrow(macromap_coloc),
    contexts = n_distinct(macromap_mr_long$context),
    unique_genes = n_distinct(macromap_mr_long$gene_symbol)
  )
)
write_csv(input_audit, file.path(dirs$tables, "input_audit_summary.csv"))

cat("\nInput parsing completed.\n")
print(input_audit)
cat("MacroMap Ensembl mapping rate: ",
    percent(mean(macro_gene_map$mapping_status == "Mapped_to_symbol")), "\n")
cat("MR representative selection mode: ", mr_method_selection_mode, "\n")


# =========================
# 5. Two-stage evidence table
# =========================
mr_all <- mr_all %>%
  group_by(source, context) %>%
  mutate(mr_fdr_context = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  mutate(mr_selected = pval < mr_screen_p)

ev_all <- mr_all %>%
  left_join(
    coloc_all %>% select(source, gene_symbol, context, PP.H3, PP.H4, lead_snp, n_coloc_snps),
    by = c("source", "gene_symbol", "context")
  ) %>%
  left_join(context_map, by = c("source", "context")) %>%
  mutate(
    coloc_success = mr_selected & is.finite(PP.H4),
    unexpected_coloc = !mr_selected & is.finite(PP.H4),
    raw_score = case_when(
      coloc_success ~ pmin(-log10(pmax(pval, min_p)), max_logp) * PP.H4,
      TRUE ~ NA_real_
    ),
    evidence_status = case_when(
      !mr_selected ~ "Not_selected_for_coloc",
      mr_selected & !is.finite(PP.H4) ~ "Selected_but_coloc_unavailable",
      PP.H4 >= strong_coloc_cutoff ~ "Strong_coloc",
      PP.H4 >= moderate_coloc_cutoff ~ "Moderate_coloc",
      TRUE ~ "Weak_coloc"
    )
  ) %>%
  group_by(source, context) %>%
  mutate(
    conditional_pct_score = rank_percentile(raw_score),
    signed_conditional_pct = sign(beta) * conditional_pct_score,
    n_coloc_rank_reference = sum(is.finite(raw_score))
  ) %>%
  ungroup() %>%
  mutate(
    # 0 means the pair did not pass the prespecified MR screening gate.
    # It does NOT mean PP.H4 = 0 or proven absence of colocalization.
    integrated_score = case_when(
      coloc_success ~ conditional_pct_score,
      !mr_selected ~ 0,
      TRUE ~ NA_real_
    ),
    signed_integrated_score = case_when(
      is.finite(integrated_score) ~ sign(beta) * integrated_score,
      TRUE ~ NA_real_
    ),
    raw_integrated_score = case_when(
      coloc_success ~ raw_score,
      !mr_selected ~ 0,
      TRUE ~ NA_real_
    ),
    signed_raw_integrated_score = case_when(
      is.finite(raw_integrated_score) ~ sign(beta) * raw_integrated_score,
      TRUE ~ NA_real_
    ),
    # Tier evidence uses nominal MR P only; mr_fdr_context is retained solely
    # as a descriptive column and never enters tier assignment.
    strong_support = coloc_success & pval < tier_mr_p & PP.H4 >= strong_coloc_cutoff,
    moderate_support = coloc_success & pval < tier_mr_p & PP.H4 >= moderate_coloc_cutoff,
    exact_context = paste(source, context, sep = "__")
  )

write_csv(ev_all, file.path(dirs$tables, "unified_two_stage_evidence.csv"))

# =========================
# 6. QC summaries and figures
# =========================
context_qc <- ev_all %>%
  group_by(source, block, module, context) %>%
  summarise(
    n_mr = n(),
    n_genes = n_distinct(gene_symbol),
    n_mr_selected = sum(mr_selected, na.rm = TRUE),
    mr_selection_rate = mean(mr_selected, na.rm = TRUE),
    n_coloc_success = sum(coloc_success, na.rm = TRUE),
    coloc_completion_among_selected = if_else(
      n_mr_selected > 0,
      n_coloc_success / n_mr_selected,
      NA_real_
    ),
    n_strong_coloc = sum(strong_support, na.rm = TRUE),
    n_moderate_coloc = sum(moderate_support, na.rm = TRUE),
    median_raw_score = median(raw_score, na.rm = TRUE),
    rank_reference_n = max(n_coloc_rank_reference, na.rm = TRUE),
    .groups = "drop"
  )

pipeline_qc <- ev_all %>%
  group_by(source) %>%
  summarise(
    MR_pairs = n(),
    MR_selected = sum(mr_selected, na.rm = TRUE),
    Coloc_success = sum(coloc_success, na.rm = TRUE),
    Moderate_coloc = sum(moderate_support, na.rm = TRUE),
    Strong_coloc = sum(strong_support, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-source, names_to = "stage", values_to = "n") %>%
  mutate(
    stage = factor(
      stage,
      levels = c("MR_pairs", "MR_selected", "Coloc_success", "Moderate_coloc", "Strong_coloc")
    )
  )

write.xlsx(
  list(Context_QC = context_qc, Pipeline_QC = pipeline_qc),
  file.path(dirs$tables, "QC_summary.xlsx"),
  overwrite = TRUE
)

p_pipeline <- ggplot(pipeline_qc, aes(stage, n, group = source, color = source)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = comma(n)), vjust = -0.6, size = 3) +
  theme_bw(base_size = 12) +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "Number of gene-context pairs", color = "Resource") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot_both(p_pipeline, file.path(dirs$qc, "01_pipeline_attrition"), 9, 6)

p_selection <- context_qc %>%
  mutate(context_label = paste(source, context, sep = ": ")) %>%
  ggplot(aes(mr_selection_rate, reorder(context_label, mr_selection_rate))) +
  geom_col() +
  facet_wrap(~source, scales = "free_y", ncol = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw(base_size = 11) +
  labs(x = "Proportion with nominal MR P < 0.05", y = NULL)
save_plot_both(p_selection, file.path(dirs$qc, "02_context_MR_selection_rate"), 10, 10)

p_completion <- context_qc %>%
  mutate(context_label = paste(source, context, sep = ": ")) %>%
  ggplot(aes(coloc_completion_among_selected, reorder(context_label, coloc_completion_among_selected))) +
  geom_col() +
  facet_wrap(~source, scales = "free_y", ncol = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(x = "Colocalization completion among MR-selected pairs", y = NULL)
save_plot_both(p_completion, file.path(dirs$qc, "03_coloc_completion_rate"), 10, 10)

p_raw_dist <- ev_all %>%
  filter(coloc_success) %>%
  ggplot(aes(source, raw_score)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = expression(-log[10](P[MR]) %*% PP.H4))
save_plot_both(p_raw_dist, file.path(dirs$qc, "04_raw_score_by_resource"), 7, 6)

p_pct_dist <- ev_all %>%
  filter(coloc_success) %>%
  ggplot(aes(conditional_pct_score, fill = source)) +
  geom_histogram(bins = 20, position = "identity", alpha = 0.45) +
  facet_wrap(~source, ncol = 1) +
  theme_bw(base_size = 12) +
  labs(x = "Conditional percentile score", y = "Count", fill = "Resource")
save_plot_both(p_pct_dist, file.path(dirs$qc, "05_percentile_score_distribution"), 8, 7)

# =========================
# 7. Candidate gene universe
# =========================
gene_availability <- ev_all %>%
  group_by(gene_symbol) %>%
  summarise(
    has_static_mr = any(source == "OneK1K_static"),
    has_dynamic_mr = any(source == "Macrophage_dynamic"),
    any_mr_selected = any(mr_selected, na.rm = TRUE),
    any_coloc_success = any(coloc_success, na.rm = TRUE),
    static_coloc_success = any(source == "OneK1K_static" & coloc_success, na.rm = TRUE),
    dynamic_coloc_success = any(source == "Macrophage_dynamic" & coloc_success, na.rm = TRUE),
    .groups = "drop"
  )

candidate_genes <- gene_availability %>%
  filter(
    case_when(
      candidate_mode == "common_any_coloc" ~ has_static_mr & has_dynamic_mr & any_coloc_success,
      candidate_mode == "common_any_selected" ~ has_static_mr & has_dynamic_mr & any_mr_selected,
      candidate_mode == "union_any_coloc" ~ any_coloc_success,
      TRUE ~ FALSE
    )
  ) %>%
  pull(gene_symbol)

if (length(candidate_genes) < 6) {
  stop("Fewer than six candidate genes remain. Review candidate_mode and input coverage.")
}

write_csv(
  gene_availability %>% mutate(in_joint_candidate_set = gene_symbol %in% candidate_genes),
  file.path(dirs$tables, "gene_availability_and_candidate_status.csv")
)

# =========================
# 8. Module aggregation
# =========================
module_reference <- context_map %>%
  group_by(block, source, module, feature) %>%
  summarise(n_expected_contexts = n_distinct(context), .groups = "drop")

aggregate_modules <- function(ev, genes, aggregation_method) {
  summary_long <- ev %>%
    filter(gene_symbol %in% genes) %>%
    group_by(gene_symbol, source, block, module, feature) %>%
    summarise(
      n_available_contexts = n_distinct(context),
      n_mr_selected = sum(mr_selected, na.rm = TRUE),
      n_coloc_success = sum(coloc_success, na.rm = TRUE),
      module_score_raw = aggregate_vector(integrated_score, aggregation_method, signed = FALSE),
      module_signed_raw = aggregate_vector(signed_integrated_score, aggregation_method, signed = TRUE),
      support_breadth_moderate = mean(moderate_support, na.rm = TRUE),
      support_breadth_strong = mean(strong_support, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(
      crossing(tibble(gene_symbol = genes), module_reference),
      by = c("gene_symbol", "source", "block", "module", "feature")
    ) %>%
    mutate(
      n_available_contexts = replace_na(n_available_contexts, 0L),
      n_mr_selected = replace_na(n_mr_selected, 0L),
      n_coloc_success = replace_na(n_coloc_success, 0L),
      context_coverage = n_available_contexts / n_expected_contexts,
      module_score = if_else(
        context_coverage >= min_context_coverage_within_module,
        module_score_raw,
        NA_real_
      ),
      module_signed_score = if_else(
        context_coverage >= min_context_coverage_within_module,
        module_signed_raw,
        NA_real_
      )
    )

  score_mat <- summary_long %>%
    select(gene_symbol, feature, module_score) %>%
    pivot_wider(names_from = feature, values_from = module_score) %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()

  signed_mat <- summary_long %>%
    select(gene_symbol, feature, module_signed_score) %>%
    pivot_wider(names_from = feature, values_from = module_signed_score) %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()

  list(long = summary_long, score = score_mat, signed = signed_mat)
}

main_modules <- aggregate_modules(ev_all, candidate_genes, main_aggregation)

# Block-level module coverage filter
module_coverage_gene <- main_modules$long %>%
  group_by(gene_symbol, block) %>%
  summarise(module_coverage = mean(is.finite(module_score)), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = module_coverage, names_prefix = "coverage_")

main_genes <- module_coverage_gene %>%
  mutate(
    coverage_Static = replace_na(coverage_Static, 0),
    coverage_Dynamic = replace_na(coverage_Dynamic, 0)
  ) %>%
  filter(
    coverage_Static >= min_module_coverage_per_block,
    coverage_Dynamic >= min_module_coverage_per_block
  ) %>%
  pull(gene_symbol)

if (length(main_genes) < 6) {
  stop(
    "Fewer than six genes pass block coverage filters. Consider reducing ",
    "min_module_coverage_per_block or reviewing context mapping."
  )
}

write_csv(module_coverage_gene, file.path(dirs$tables, "gene_block_module_coverage.csv"))
write_csv(main_modules$long, file.path(dirs$tables, "module_aggregation_main.csv"))

# =========================
# 9. Block standardization and weighting
# =========================
prepare_block_matrix <- function(module_bundle, genes, static_weight, dynamic_weight) {
  mat <- module_bundle$score[genes, , drop = FALSE]
  signed_mat <- module_bundle$signed[genes, , drop = FALSE]

  static_cols <- grep("^Static__", colnames(mat), value = TRUE)
  dynamic_cols <- grep("^Dynamic__", colnames(mat), value = TRUE)

  if (length(static_cols) == 0 || length(dynamic_cols) == 0) {
    stop("Both Static and Dynamic modules are required for joint clustering.")
  }

  static_raw <- mat[, static_cols, drop = FALSE]
  dynamic_raw <- mat[, dynamic_cols, drop = FALSE]

  keep_static <- apply(static_raw, 2, sd, na.rm = TRUE) > 0
  keep_dynamic <- apply(dynamic_raw, 2, sd, na.rm = TRUE) > 0
  keep_static[is.na(keep_static)] <- FALSE
  keep_dynamic[is.na(keep_dynamic)] <- FALSE

  static_raw <- static_raw[, keep_static, drop = FALSE]
  dynamic_raw <- dynamic_raw[, keep_dynamic, drop = FALSE]

  if (ncol(static_raw) == 0 || ncol(dynamic_raw) == 0) {
    stop("No variable module remains in one of the two data blocks after filtering.")
  }

  static_z <- scale(static_raw)
  dynamic_z <- scale(dynamic_raw)
  static_z[is.na(static_z)] <- 0
  dynamic_z[is.na(dynamic_z)] <- 0

  static_balanced <- static_z * sqrt(static_weight / ncol(static_z))
  dynamic_balanced <- dynamic_z * sqrt(dynamic_weight / ncol(dynamic_z))

  joint <- cbind(static_balanced, dynamic_balanced)
  profile <- cbind(static_z, dynamic_z)

  signed_selected <- signed_mat[, c(colnames(static_raw), colnames(dynamic_raw)), drop = FALSE]
  signed_selected[is.na(signed_selected)] <- 0

  list(
    joint = joint,
    profile = profile,
    signed = signed_selected,
    static_z = static_z,
    dynamic_z = dynamic_z,
    static_raw = static_raw,
    dynamic_raw = dynamic_raw,
    static_weight = static_weight,
    dynamic_weight = dynamic_weight
  )
}

main_blocks <- prepare_block_matrix(
  main_modules,
  main_genes,
  main_static_weight,
  main_dynamic_weight
)

write.xlsx(
  list(
    Joint_weighted = as.data.frame(main_blocks$joint) %>% rownames_to_column("gene_symbol"),
    Profile_unweighted_Z = as.data.frame(main_blocks$profile) %>% rownames_to_column("gene_symbol"),
    Signed_module = as.data.frame(main_blocks$signed) %>% rownames_to_column("gene_symbol"),
    Static_raw = as.data.frame(main_blocks$static_raw) %>% rownames_to_column("gene_symbol"),
    Dynamic_raw = as.data.frame(main_blocks$dynamic_raw) %>% rownames_to_column("gene_symbol")
  ),
  file.path(dirs$tables, "main_joint_matrices.xlsx"),
  overwrite = TRUE
)

# =========================
# 10. Consensus clustering function
# =========================
run_consensus <- function(mat, analysis_name, out_path, distance_method, reps) {
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  n_gene <- nrow(mat)
  if (n_gene < 6) stop("Too few genes for clustering: ", analysis_name)

  maxK <- min(max_k_cap, max(2, floor(n_gene / 2)))
  linkage <- ifelse(distance_method == "euclidean", "ward.D2", "average")

  set.seed(12345)
  ccp <- ConsensusClusterPlus(
    d = t(mat),
    maxK = maxK,
    reps = reps,
    pItem = p_item,
    pFeature = min(p_feature, 1),
    clusterAlg = "hc",
    distance = distance_method,
    innerLinkage = linkage,
    finalLinkage = linkage,
    seed = 12345,
    plot = "png",
    title = file.path(out_path, "ConsensusClusterPlus")
  )

  pac_tbl <- tibble(
    K = 2:maxK,
    PAC = map_dbl(2:maxK, ~calc_pac(ccp[[.x]]$consensusMatrix))
  )

  cluster_size_tbl <- map_dfr(2:maxK, function(k) {
    tibble(
      K = k,
      gene_symbol = names(ccp[[k]]$consensusClass),
      cluster = as.integer(ccp[[k]]$consensusClass)
    ) %>%
      count(K, cluster, name = "n")
  })

  valid_k <- cluster_size_tbl %>%
    group_by(K) %>%
    summarise(min_cluster_n = min(n), .groups = "drop") %>%
    filter(min_cluster_n >= min_cluster_size)

  if (nrow(valid_k) == 0) {
    bestK <- 2
  } else {
    bestK <- pac_tbl %>%
      filter(K %in% valid_k$K) %>%
      arrange(PAC, K) %>%
      slice(1) %>%
      pull(K)
  }

  cluster_res <- tibble(
    gene_symbol = names(ccp[[bestK]]$consensusClass),
    cluster_id = as.integer(ccp[[bestK]]$consensusClass),
    cluster = paste0("Cluster_", cluster_id)
  )

  dist_obj <- make_dist(mat[cluster_res$gene_symbol, , drop = FALSE], distance_method)
  sil <- silhouette(cluster_res$cluster_id, dist_obj)
  sil_df <- as.data.frame(sil) %>%
    rownames_to_column("gene_symbol") %>%
    rename(
      cluster_id = cluster,
      neighbor_cluster = neighbor,
      silhouette_width = sil_width
    )

  p_pac <- ggplot(pac_tbl, aes(K, PAC)) +
    geom_line() + geom_point(size = 2.5) +
    geom_vline(xintercept = bestK, linetype = 2) +
    theme_bw(base_size = 12) +
    scale_x_continuous(breaks = 2:maxK) +
    labs(x = "Number of clusters (K)", y = "PAC")
  save_plot_both(p_pac, file.path(out_path, "PAC_by_K"), 7, 5)

  png(file.path(out_path, "silhouette_plot.png"), width = 1800, height = 1300, res = 220)
  plot(sil, border = NA, main = paste0(analysis_name, "; K = ", bestK))
  dev.off()
  pdf(file.path(out_path, "silhouette_plot.pdf"), width = 9, height = 6.5)
  plot(sil, border = NA, main = paste0(analysis_name, "; K = ", bestK))
  dev.off()

  pca <- prcomp(mat[cluster_res$gene_symbol, , drop = FALSE], center = FALSE, scale. = FALSE)
  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("gene_symbol") %>%
    left_join(cluster_res, by = "gene_symbol")

  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cluster, label = gene_symbol)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 60) +
    theme_bw(base_size = 12) +
    labs(color = "Cluster")
  save_plot_both(p_pca, file.path(out_path, "PCA_clusters"), 9, 7)

  write.xlsx(
    list(
      Cluster_assignment = cluster_res,
      PAC = pac_tbl,
      Cluster_size = cluster_size_tbl,
      Silhouette = sil_df,
      PCA = pca_df
    ),
    file.path(out_path, "clustering_results.xlsx"),
    overwrite = TRUE
  )

  list(
    name = analysis_name,
    ccp = ccp,
    bestK = bestK,
    clusters = cluster_res,
    pac = pac_tbl,
    cluster_size = cluster_size_tbl,
    silhouette = sil_df,
    pca = pca_df,
    mat = mat,
    distance = distance_method
  )
}

# =========================
# 11. Main joint consensus clustering
# =========================
main_cluster <- run_consensus(
  main_blocks$joint,
  analysis_name = "DynaMac-GP block-balanced joint clustering",
  out_path = dirs$main,
  distance_method = main_distance,
  reps = main_reps
)

cluster_res <- main_cluster$clusters

# =========================
# 12. Main heatmaps and cluster profiles
# =========================
feature_annotation <- tibble(feature = colnames(main_blocks$profile)) %>%
  separate(feature, into = c("block", "module"), sep = "__", remove = FALSE) %>%
  column_to_rownames("feature")

cluster_order <- cluster_res %>% arrange(cluster_id, gene_symbol)
ordered_genes <- cluster_order$gene_symbol
cluster_sizes <- table(cluster_order$cluster_id)
gaps_row <- cumsum(as.integer(cluster_sizes))
gaps_row <- gaps_row[-length(gaps_row)]

ann_row <- cluster_order %>%
  select(gene_symbol, cluster) %>%
  column_to_rownames("gene_symbol")

plot_heatmap <- function(mat, filename, main_text, width = 12, height = 10) {
  png(paste0(filename, ".png"), width = 2800, height = 2300, res = 220)
  pheatmap(
    mat[ordered_genes, , drop = FALSE],
    annotation_row = ann_row,
    annotation_col = feature_annotation[colnames(mat), , drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    gaps_row = gaps_row,
    show_rownames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    main = main_text
  )
  dev.off()

  pdf(paste0(filename, ".pdf"), width = width, height = height)
  pheatmap(
    mat[ordered_genes, , drop = FALSE],
    annotation_row = ann_row,
    annotation_col = feature_annotation[colnames(mat), , drop = FALSE],
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    gaps_row = gaps_row,
    show_rownames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    main = main_text
  )
  dev.off()
}

plot_heatmap(
  main_blocks$joint,
  file.path(dirs$main, "main_joint_weighted_heatmap"),
  paste0("Block-balanced joint evidence; K = ", main_cluster$bestK)
)

plot_heatmap(
  main_blocks$profile,
  file.path(dirs$main, "main_unweighted_profile_heatmap"),
  "Unweighted standardized static-dynamic profiles"
)

plot_heatmap(
  main_blocks$signed,
  file.path(dirs$main, "main_signed_module_heatmap"),
  "Signed conditional-percentile evidence"
)

cluster_profile <- as.data.frame(main_blocks$profile) %>%
  rownames_to_column("gene_symbol") %>%
  left_join(cluster_res, by = "gene_symbol") %>%
  pivot_longer(
    cols = -c(gene_symbol, cluster_id, cluster),
    names_to = "feature",
    values_to = "z_score"
  ) %>%
  separate(feature, into = c("block", "module"), sep = "__", remove = FALSE) %>%
  group_by(cluster_id, cluster, block, module, feature) %>%
  summarise(
    mean_z = mean(z_score),
    median_z = median(z_score),
    mean_abs_z = mean(abs(z_score)),
    .groups = "drop"
  )

centroid_mat <- cluster_profile %>%
  select(cluster, feature, mean_z) %>%
  pivot_wider(names_from = feature, values_from = mean_z) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

png(file.path(dirs$main, "cluster_centroid_heatmap.png"), width = 2400, height = 1400, res = 220)
pheatmap(
  centroid_mat,
  annotation_col = feature_annotation[colnames(centroid_mat), , drop = FALSE],
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize_col = 8,
  main = "Cluster centroids across static and dynamic modules"
)
dev.off()
pdf(file.path(dirs$main, "cluster_centroid_heatmap.pdf"), width = 12, height = 7)
pheatmap(
  centroid_mat,
  annotation_col = feature_annotation[colnames(centroid_mat), , drop = FALSE],
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize_col = 8,
  main = "Cluster centroids across static and dynamic modules"
)
dev.off()

cluster_labels <- cluster_profile %>%
  group_by(cluster) %>%
  summarise(
    top_static = feature[block == "Static"][which.max(mean_z[block == "Static"])],
    top_static_z = max(mean_z[block == "Static"], na.rm = TRUE),
    top_dynamic = feature[block == "Dynamic"][which.max(mean_z[block == "Dynamic"])],
    top_dynamic_z = max(mean_z[block == "Dynamic"], na.rm = TRUE),
    static_contribution = mean(mean_abs_z[block == "Static"], na.rm = TRUE),
    dynamic_contribution = mean(mean_abs_z[block == "Dynamic"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    top_static = str_remove(top_static, "^Static__"),
    top_dynamic = str_remove(top_dynamic, "^Dynamic__"),
    automatic_label = case_when(
      top_static_z > 0.25 & top_dynamic_z > 0.25 ~
        paste0("Cross-resource: ", top_static, " + ", top_dynamic),
      top_static_z > 0.25 & top_dynamic_z <= 0.25 ~
        paste0("Static-dominant: ", top_static),
      top_dynamic_z > 0.25 & top_static_z <= 0.25 ~
        paste0("Dynamic-dominant: ", top_dynamic),
      TRUE ~ "Low/neutral joint-evidence pattern"
    )
  )

block_contribution_long <- cluster_labels %>%
  select(cluster, static_contribution, dynamic_contribution) %>%
  pivot_longer(-cluster, names_to = "block", values_to = "mean_absolute_z") %>%
  mutate(block = recode(
    block,
    static_contribution = "Static block",
    dynamic_contribution = "Dynamic block"
  ))

p_block <- ggplot(block_contribution_long, aes(cluster, mean_absolute_z, fill = block)) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Mean absolute standardized contribution", fill = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot_both(p_block, file.path(dirs$main, "cluster_block_contributions"), 9, 6)

write.xlsx(
  list(
    Cluster_profile = cluster_profile,
    Automatic_labels = cluster_labels,
    Centroid_matrix = as.data.frame(centroid_mat) %>% rownames_to_column("cluster")
  ),
  file.path(dirs$main, "cluster_profiles_and_labels.xlsx"),
  overwrite = TRUE
)

# =========================
# 13. Source-specific nominal-P tiering
# =========================
best_source_row <- function(df, source_value) {
  df %>%
    filter(source == source_value, coloc_success) %>%
    group_by(gene_symbol) %>%
    arrange(desc(raw_score), desc(PP.H4), pval, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene_symbol,
      best_context = context,
      best_module = module,
      best_beta = beta,
      best_pval = pval,
      best_fdr = mr_fdr_context,
      best_PP.H4 = PP.H4,
      best_raw_score = raw_score,
      best_pct_score = conditional_pct_score
    )
}

static_best <- best_source_row(ev_all, "OneK1K_static") %>%
  rename_with(~paste0("static_", .x), -gene_symbol)

dynamic_best <- best_source_row(ev_all, "Macrophage_dynamic") %>%
  rename_with(~paste0("dynamic_", .x), -gene_symbol)

gene_support <- ev_all %>%
  group_by(gene_symbol) %>%
  summarise(
    n_static_mr = sum(source == "OneK1K_static"),
    n_dynamic_mr = sum(source == "Macrophage_dynamic"),
    n_static_selected = sum(source == "OneK1K_static" & mr_selected, na.rm = TRUE),
    n_dynamic_selected = sum(source == "Macrophage_dynamic" & mr_selected, na.rm = TRUE),
    n_static_coloc = sum(source == "OneK1K_static" & coloc_success, na.rm = TRUE),
    n_dynamic_coloc = sum(source == "Macrophage_dynamic" & coloc_success, na.rm = TRUE),

    static_strong = any(
      source == "OneK1K_static" & strong_support,
      na.rm = TRUE
    ),
    static_moderate = any(
      source == "OneK1K_static" & moderate_support,
      na.rm = TRUE
    ),
    dynamic_strong = any(
      source == "Macrophage_dynamic" & strong_support,
      na.rm = TRUE
    ),
    dynamic_moderate = any(
      source == "Macrophage_dynamic" & moderate_support,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(static_best, by = "gene_symbol") %>%
  left_join(dynamic_best, by = "gene_symbol") %>%
  mutate(
    # Final tier system:
    # Tier 1: both resources show strong evidence.
    # Tier 2A: dynamic strong + static moderate.
    # Tier 2B: static strong + dynamic moderate.
    # Tier 3A: dynamic-only strong evidence.
    # Tier 3B: static-only strong evidence.
    # Tier 4: both resources show moderate evidence.
    # Tier 5: all remaining patterns.
    tier = case_when(
      static_strong & dynamic_strong ~
        "Tier1_cross_resource_strong",
      dynamic_strong & static_moderate ~
        "Tier2A_dynamic_strong_with_static_support",
      static_strong & dynamic_moderate ~
        "Tier2B_static_strong_with_dynamic_support",
      dynamic_strong ~
        "Tier3A_dynamic_specific_strong",
      static_strong ~
        "Tier3B_static_specific_strong",
      static_moderate & dynamic_moderate ~
        "Tier4_cross_resource_moderate",
      TRUE ~
        "Tier5_limited_evidence"
    ),
    tier_rank = case_when(
      str_detect(tier, "^Tier1") ~ 1L,
      str_detect(tier, "^Tier2A") ~ 2L,
      str_detect(tier, "^Tier2B") ~ 3L,
      str_detect(tier, "^Tier3A") ~ 4L,
      str_detect(tier, "^Tier3B") ~ 5L,
      str_detect(tier, "^Tier4") ~ 6L,
      TRUE ~ 7L
    ),
    evidence_class = case_when(
      str_detect(tier, "^Tier1|^Tier2|^Tier4") ~
        "Cross-resource convergent evidence",
      str_detect(tier, "^Tier3A") ~
        "Dynamic-specific strong evidence",
      str_detect(tier, "^Tier3B") ~
        "Static-specific strong evidence",
      TRUE ~
        "Limited evidence"
    ),
    direction_concordant = case_when(
      is.finite(static_best_beta) & is.finite(dynamic_best_beta) ~
        sign(static_best_beta) == sign(dynamic_best_beta),
      TRUE ~ NA
    ),
    integrated_rank_score =
      0.5 * replace_na(static_best_pct_score, 0) +
      0.5 * replace_na(dynamic_best_pct_score, 0)
  ) %>%
  left_join(cluster_res, by = "gene_symbol") %>%
  left_join(
    cluster_labels %>% select(cluster, automatic_label),
    by = "cluster"
  ) %>%
  arrange(tier_rank, desc(integrated_rank_score))

tier_definitions <- tribble(
  ~tier, ~tier_rank, ~criteria, ~interpretation,
  "Tier1_cross_resource_strong", 1L,
  "Static strong AND dynamic strong",
  "Strong MR-colocalization evidence in both OneK1K and MacroMap",
  "Tier2A_dynamic_strong_with_static_support", 2L,
  "Dynamic strong AND static moderate",
  "Dynamic macrophage evidence is strong, with supporting static immune evidence",
  "Tier2B_static_strong_with_dynamic_support", 3L,
  "Static strong AND dynamic moderate",
  "Static immune evidence is strong, with supporting dynamic macrophage evidence",
  "Tier3A_dynamic_specific_strong", 4L,
  "Dynamic strong; static does not reach moderate evidence",
  "Dynamic macrophage-specific strong evidence",
  "Tier3B_static_specific_strong", 5L,
  "Static strong; dynamic does not reach moderate evidence",
  "Static immune-cell-specific strong evidence",
  "Tier4_cross_resource_moderate", 6L,
  "Static moderate AND dynamic moderate",
  "Moderate convergent evidence across both resources",
  "Tier5_limited_evidence", 7L,
  "All other evidence combinations",
  "Limited, weak, or single-resource moderate evidence"
)

write_csv(
  tier_definitions,
  file.path(dirs$tier, "tier_definitions.csv")
)

write_csv(gene_support, file.path(dirs$tier, "gene_tiers_and_clusters.csv"))

p_tier <- gene_support %>%
  count(tier, tier_rank, name = "n") %>%
  ggplot(aes(n, reorder(tier, -tier_rank))) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  theme_bw(base_size = 11) +
  expand_limits(x = max(table(gene_support$tier)) * 1.12) +
  labs(x = "Number of genes", y = NULL)
save_plot_both(p_tier, file.path(dirs$tier, "01_tier_counts"), 10, 6.5)

p_cluster_tier <- gene_support %>%
  filter(!is.na(cluster)) %>%
  count(cluster, tier) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(cluster, prop, fill = tier)) +
  geom_col() +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Proportion within cluster", fill = "Tier") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot_both(p_cluster_tier, file.path(dirs$tier, "02_cluster_by_tier"), 10, 7)

direction_summary <- gene_support %>%
  filter(
    static_moderate,
    dynamic_moderate,
    !is.na(direction_concordant)
  ) %>%
  count(direction_concordant, name = "n") %>%
  mutate(
    direction = if_else(direction_concordant, "Concordant", "Discordant"),
    proportion = n / sum(n)
  )

if (nrow(direction_summary) > 0) {
  p_direction <- ggplot(direction_summary, aes(direction, proportion, fill = direction)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = percent(proportion, accuracy = 0.1)), vjust = -0.4) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1.08)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(x = NULL, y = "Proportion among cross-resource supported genes")
  save_plot_both(p_direction, file.path(dirs$tier, "03_direction_concordance"), 7, 5.5)
}

write.xlsx(
  list(
    Tier_definitions = tier_definitions,
    Gene_tiers = gene_support,
    Direction_summary = direction_summary,
    Tier_cluster_crosstab = gene_support %>% count(cluster, tier)
  ),
  file.path(dirs$tier, "tier_summary.xlsx"),
  overwrite = TRUE
)

# =========================
# 14. Sensitivity analyses
# =========================
sensitivity_configs <- tribble(
  ~analysis, ~aggregation, ~static_weight, ~dynamic_weight, ~distance, ~block_mode,
  "Main_50S_50D_topQ_euclidean", "top_quartile", 0.50, 0.50, "euclidean", "joint",
  "Weight_30S_70D",              "top_quartile", 0.30, 0.70, "euclidean", "joint",
  "Weight_70S_30D",              "top_quartile", 0.70, 0.30, "euclidean", "joint",
  "Aggregation_max",             "max",          0.50, 0.50, "euclidean", "joint",
  "Aggregation_median",          "median",       0.50, 0.50, "euclidean", "joint",
  "Distance_pearson",            "top_quartile", 0.50, 0.50, "pearson",   "joint",
  "Static_only",                 "top_quartile", 1.00, 0.00, "euclidean", "static",
  "Dynamic_only",                "top_quartile", 0.00, 1.00, "euclidean", "dynamic"
)

sensitivity_results <- list()

for (i in seq_len(nrow(sensitivity_configs))) {
  cfg <- sensitivity_configs[i, ]
  analysis_name <- cfg$analysis
  message("Running sensitivity analysis: ", analysis_name)

  if (analysis_name == "Main_50S_50D_topQ_euclidean") {
    sensitivity_results[[analysis_name]] <- main_cluster
    next
  }

  bundle <- if (cfg$aggregation == main_aggregation) {
    main_modules
  } else {
    aggregate_modules(ev_all, candidate_genes, cfg$aggregation)
  }

  blocks <- prepare_block_matrix(
    bundle,
    main_genes,
    static_weight = ifelse(cfg$block_mode == "dynamic", 0.0001, max(cfg$static_weight, 0.0001)),
    dynamic_weight = ifelse(cfg$block_mode == "static", 0.0001, max(cfg$dynamic_weight, 0.0001))
  )

  sens_mat <- if (cfg$block_mode == "static") {
    blocks$static_z
  } else if (cfg$block_mode == "dynamic") {
    blocks$dynamic_z
  } else {
    blocks$joint
  }

  result <- tryCatch(
    run_consensus(
      sens_mat,
      analysis_name = analysis_name,
      out_path = file.path(dirs$sensitivity, analysis_name),
      distance_method = cfg$distance,
      reps = sensitivity_reps
    ),
    error = function(e) {
      warning("Sensitivity analysis failed: ", analysis_name, "; ", conditionMessage(e))
      NULL
    }
  )

  sensitivity_results[[analysis_name]] <- result
}

# Diagnostic reproduction of the original-style, unharmonized score clustering.
# This is reported only to quantify how much the revised harmonization changes
# the original solution; it is not treated as the preferred analysis.
ev_raw_diagnostic <- ev_all %>%
  mutate(
    integrated_score = raw_integrated_score,
    signed_integrated_score = signed_raw_integrated_score
  )

raw_modules_diagnostic <- aggregate_modules(
  ev_raw_diagnostic,
  candidate_genes,
  aggregation_method = "max"
)

raw_diag_mat <- raw_modules_diagnostic$score[main_genes, , drop = FALSE]
raw_keep <- apply(raw_diag_mat, 2, sd, na.rm = TRUE) > 0
raw_keep[is.na(raw_keep)] <- FALSE
raw_diag_mat <- raw_diag_mat[, raw_keep, drop = FALSE]
raw_diag_scaled <- scale(raw_diag_mat)
raw_diag_scaled[is.na(raw_diag_scaled)] <- 0

raw_diagnostic_result <- tryCatch(
  run_consensus(
    raw_diag_scaled,
    analysis_name = "Legacy-style unharmonized raw-score clustering",
    out_path = file.path(dirs$sensitivity, "Legacy_raw_unharmonized"),
    distance_method = "pearson",
    reps = sensitivity_reps
  ),
  error = function(e) {
    warning("Legacy raw-score diagnostic failed: ", conditionMessage(e))
    NULL
  }
)
sensitivity_results[["Legacy_raw_unharmonized"]] <- raw_diagnostic_result

main_assignment <- main_cluster$clusters %>%
  select(gene_symbol, main_cluster = cluster)

sensitivity_comparison <- map_dfr(names(sensitivity_results), function(nm) {
  res <- sensitivity_results[[nm]]
  if (is.null(res)) {
    return(tibble(
      analysis = nm, n_common = NA_integer_, K = NA_integer_,
      ARI = NA_real_, NMI = NA_real_, coclustering_Jaccard = NA_real_
    ))
  }

  comp <- main_assignment %>%
    inner_join(
      res$clusters %>% select(gene_symbol, sensitivity_cluster = cluster),
      by = "gene_symbol"
    )

  tibble(
    analysis = nm,
    n_common = nrow(comp),
    K = res$bestK,
    ARI = adjusted_rand_index(comp$main_cluster, comp$sensitivity_cluster),
    NMI = normalized_mutual_information(comp$main_cluster, comp$sensitivity_cluster),
    coclustering_Jaccard = coclustering_jaccard(comp$main_cluster, comp$sensitivity_cluster)
  )
})

write_csv(sensitivity_comparison, file.path(dirs$sensitivity, "sensitivity_cluster_stability.csv"))

p_sens <- sensitivity_comparison %>%
  filter(analysis != "Main_50S_50D_topQ_euclidean") %>%
  select(analysis, ARI, NMI, coclustering_Jaccard) %>%
  pivot_longer(-analysis, names_to = "metric", values_to = "value") %>%
  ggplot(aes(value, reorder(analysis, value), fill = metric)) +
  geom_col(position = "dodge") +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(x = "Agreement with the main clustering", y = NULL, fill = "Metric")
save_plot_both(p_sens, file.path(dirs$sensitivity, "cluster_stability_metrics"), 10, 7)


# =========================
# 15. Optional enrichment analysis
# =========================
run_group_enrichment <- function(group_df, group_col, universe_symbols, prefix) {
  needed <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
  if (!all(map_lgl(needed, requireNamespace, quietly = TRUE))) {
    warning("Enrichment packages unavailable; enrichment skipped.")
    return(NULL)
  }

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
  })

  universe_map <- bitr(
    unique(universe_symbols),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  universe_entrez <- unique(universe_map$ENTREZID)

  groups <- split(group_df$gene_symbol, group_df[[group_col]])
  go_tables <- list()
  kegg_tables <- list()

  for (grp in names(groups)) {
    symbols <- unique(groups[[grp]])
    if (length(symbols) < 5) next

    mapped <- bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    entrez <- unique(mapped$ENTREZID)
    if (length(entrez) < 5) next

    ego <- tryCatch(
      enrichGO(
        gene = entrez,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20,
        readable = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      go_tables[[grp]] <- as.data.frame(ego)
      p <- enrichplot::dotplot(ego, showCategory = 15) + ggtitle(paste(prefix, grp, "GO-BP"))
      save_plot_both(p, file.path(dirs$enrichment, paste0(prefix, "_", make.names(grp), "_GO_BP")), 10, 7)
    }

    ekegg <- tryCatch(
      enrichKEGG(
        gene = entrez,
        universe = universe_entrez,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.20
      ),
      error = function(e) NULL
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      kegg_tables[[grp]] <- as.data.frame(ekegg)
      p <- enrichplot::dotplot(ekegg, showCategory = 15) + ggtitle(paste(prefix, grp, "KEGG"))
      save_plot_both(p, file.path(dirs$enrichment, paste0(prefix, "_", make.names(grp), "_KEGG")), 10, 7)
    }
  }

  if (length(go_tables) > 0) {
    write.xlsx(go_tables, file.path(dirs$enrichment, paste0(prefix, "_GO_BP.xlsx")), overwrite = TRUE)
  }
  if (length(kegg_tables) > 0) {
    write.xlsx(kegg_tables, file.path(dirs$enrichment, paste0(prefix, "_KEGG.xlsx")), overwrite = TRUE)
  }

  list(GO = go_tables, KEGG = kegg_tables)
}

if (run_enrichment) {
  try(
    run_group_enrichment(
      gene_support %>% filter(!is.na(cluster)),
      "cluster",
      universe_symbols = main_genes,
      prefix = "Cluster"
    ),
    silent = TRUE
  )

  tier_enrich_df <- gene_support %>%
    filter(str_detect(tier, "^Tier1|^Tier2")) %>%
    mutate(tier_group = str_extract(tier, "^Tier[0-9][A-Z]?"))

  try(
    run_group_enrichment(
      tier_enrich_df,
      "tier_group",
      universe_symbols = unique(ev_all$gene_symbol),
      prefix = "Tier"
    ),
    silent = TRUE
  )
}

# =========================
# 16. Final workbooks and RDS
# =========================
final_summary <- tibble(
  item = c(
    "MR_screen_p", "Tier_MR_p", "Strong_coloc_cutoff", "Moderate_coloc_cutoff",
    "Tier_uses_FDR", "Candidate_mode", "Main_aggregation", "Top_fraction",
    "Static_weight", "Dynamic_weight", "Main_distance",
    "Main_consensus_reps", "Candidate_genes_before_coverage",
    "Genes_in_main_clustering", "Features_in_joint_matrix", "Best_K"
  ),
  value = c(
    mr_screen_p, tier_mr_p, strong_coloc_cutoff, moderate_coloc_cutoff,
    FALSE, candidate_mode, main_aggregation, top_fraction,
    main_static_weight, main_dynamic_weight, main_distance,
    main_reps, length(candidate_genes), length(main_genes),
    ncol(main_blocks$joint), main_cluster$bestK
  )
)

write.xlsx(
  list(
    Run_summary = final_summary,
    Context_QC = context_qc,
    Candidate_status = gene_availability,
    Module_coverage = module_coverage_gene,
    Cluster_assignment = cluster_res,
    Cluster_labels = cluster_labels,
    Tier_definitions = tier_definitions,
    Gene_tiers = gene_support,
    Sensitivity = sensitivity_comparison
  ),
  file.path(out_dir, "DynaMac_GP_v2.3_final_summary.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    parameters = list(
      mr_screen_p = mr_screen_p,
      tier_mr_p = tier_mr_p,
      tier_uses_fdr = FALSE,
      strong_coloc_cutoff = strong_coloc_cutoff,
      moderate_coloc_cutoff = moderate_coloc_cutoff,
      candidate_mode = candidate_mode,
      main_aggregation = main_aggregation,
      static_weight = main_static_weight,
      dynamic_weight = main_dynamic_weight,
      distance = main_distance
    ),
    evidence_table = ev_all,
    context_qc = context_qc,
    candidate_genes = candidate_genes,
    main_genes = main_genes,
    main_modules = main_modules,
    main_blocks = main_blocks,
    main_clustering = main_cluster,
    cluster_profiles = cluster_profile,
    cluster_labels = cluster_labels,
    tier_definitions = tier_definitions,
    gene_support = gene_support,
    sensitivity_results = sensitivity_results,
    sensitivity_comparison = sensitivity_comparison
  ),
  file.path(out_dir, "DynaMac_GP_v2.3_complete_results.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("DynaMac-GP v2.3 completed\n")
cat("============================================================\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Candidate genes before coverage filter: ", length(candidate_genes), "\n", sep = "")
cat("Genes in main joint clustering: ", length(main_genes), "\n", sep = "")
cat("Static modules: ", ncol(main_blocks$static_z), "\n", sep = "")
cat("Dynamic modules: ", ncol(main_blocks$dynamic_z), "\n", sep = "")
cat("Best K: ", main_cluster$bestK, "\n", sep = "")
cat("============================================================\n")







##########FIGURE2.0#
# ============================================================
# Additional publication-quality QC figures
# Based on DynaMac_GP_v2.3_final_nominal_PP4_0.70.R
# ============================================================

library(tidyverse)
library(scales)
library(ggrepel)

# ------------------------------------------------------------
# 0. Clean labels and plotting variables
# ------------------------------------------------------------

source_labels <- c(
  "OneK1K_static" = "OneK1K",
  "Macrophage_dynamic" = "MacroMap"
)

context_qc_plot <- context_qc %>%
  mutate(
    Resource = recode(source, !!!source_labels),
    context_label = context,
    module_label = as.character(module)
  )

ev_plot <- ev_all %>%
  mutate(
    Resource = recode(source, !!!source_labels),
    log10P = pmin(
      -log10(pmax(pval, 1e-300)),
      50
    ),
    coloc_class = case_when(
      PP.H4 >= strong_coloc_cutoff ~ "Strong",
      PP.H4 >= moderate_coloc_cutoff ~ "Moderate",
      is.finite(PP.H4) ~ "Weak",
      TRUE ~ "Unavailable"
    )
  )


# ============================================================
# 1. COLOR version of original Figure 02
# Context-specific MR selection rate
# ============================================================

p_selection_color <- context_qc_plot %>%
  ggplot(
    aes(
      x = mr_selection_rate,
      y = reorder(context_label, mr_selection_rate),
      fill = module_label
    )
  ) +
  geom_col(width = 0.72) +
  facet_wrap(
    ~Resource,
    scales = "free_y",
    ncol = 1
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Gene-context pairs passing nominal MR screening (%)",
    y = NULL,
    fill = "Immune-context module"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

save_plot_both(
  p_selection_color,
  file.path(dirs$qc, "02_context_MR_selection_rate_COLOR"),
  10,
  10
)


# ============================================================
# 2. COLOR version of original Figure 03
# Colocalization completion rate
# ============================================================

p_completion_color <- context_qc_plot %>%
  filter(n_mr_selected > 0) %>%
  ggplot(
    aes(
      x = coloc_completion_among_selected,
      y = reorder(context_label, coloc_completion_among_selected),
      fill = module_label
    )
  ) +
  geom_col(width = 0.72) +
  facet_wrap(
    ~Resource,
    scales = "free_y",
    ncol = 1
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Colocalization completion among MR-screened pairs (%)",
    y = NULL,
    fill = "Immune-context module"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

save_plot_both(
  p_completion_color,
  file.path(dirs$qc, "03_coloc_completion_rate_COLOR"),
  10,
  10
)


# ============================================================
# 3. Number of moderate and strong coloc signals by context
# ============================================================

context_support_long <- context_qc_plot %>%
  select(
    Resource,
    context_label,
    module_label,
    n_moderate_coloc,
    n_strong_coloc
  ) %>%
  pivot_longer(
    cols = c(n_moderate_coloc, n_strong_coloc),
    names_to = "Evidence",
    values_to = "N"
  ) %>%
  mutate(
    Evidence = recode(
      Evidence,
      n_moderate_coloc = paste0("PP.H4 ≥ ", moderate_coloc_cutoff),
      n_strong_coloc = paste0("PP.H4 ≥ ", strong_coloc_cutoff)
    )
  )

p_support_counts <- ggplot(
  context_support_long,
  aes(
    x = N,
    y = reorder(context_label, N),
    fill = Evidence
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72
  ) +
  facet_wrap(
    ~Resource,
    scales = "free_y",
    ncol = 1
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Number of gene-context pairs",
    y = NULL,
    fill = "Colocalization support"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

save_plot_both(
  p_support_counts,
  file.path(dirs$qc, "06_context_coloc_support_counts"),
  10,
  10
)


# ============================================================
# 4. MR screening rate versus coloc completion rate
# Useful QC: does selection rate affect downstream completion?
# ============================================================

p_screen_vs_coloc <- context_qc_plot %>%
  filter(
    n_mr_selected > 0,
    is.finite(mr_selection_rate),
    is.finite(coloc_completion_among_selected)
  ) %>%
  ggplot(
    aes(
      x = mr_selection_rate,
      y = coloc_completion_among_selected,
      size = n_mr_selected,
      color = Resource
    )
  ) +
  geom_point(alpha = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = context_label),
    size = 3,
    max.overlaps = 50
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_size_continuous(
    range = c(2.5, 8)
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "MR screening rate",
    y = "Colocalization completion rate",
    color = "Resource",
    size = "MR-screened\npairs"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

save_plot_both(
  p_screen_vs_coloc,
  file.path(dirs$qc, "07_MR_selection_vs_coloc_completion"),
  9,
  7
)


# ============================================================
# 5. COLOR raw-score distribution by resource
# Replaces original black/grey violin plot
# ============================================================

p_raw_dist_color <- ev_plot %>%
  filter(
    coloc_success,
    is.finite(raw_score)
  ) %>%
  ggplot(
    aes(
      x = Resource,
      y = raw_score,
      fill = Resource
    )
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.65,
    width = 0.8
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = expression(
      RawScore == min(-log[10](P[MR]), 50) %*% PP.H4
    ),
    fill = "Resource"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

save_plot_both(
  p_raw_dist_color,
  file.path(dirs$qc, "04_raw_score_by_resource_COLOR"),
  7,
  6
)


# ============================================================
# 6. PP.H4 distribution by resource
# ============================================================

p_pph4 <- ev_plot %>%
  filter(
    coloc_success,
    is.finite(PP.H4)
  ) %>%
  ggplot(
    aes(
      x = PP.H4,
      fill = Resource
    )
  ) +
  geom_histogram(
    bins = 30,
    alpha = 0.55,
    position = "identity"
  ) +
  geom_vline(
    xintercept = c(
      moderate_coloc_cutoff,
      strong_coloc_cutoff
    ),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~Resource,
    ncol = 1,
    scales = "free_y"
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Posterior probability of shared signal (PP.H4)",
    y = "Number of gene-context pairs",
    fill = "Resource"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

save_plot_both(
  p_pph4,
  file.path(dirs$qc, "08_PPH4_distribution_by_resource"),
  8,
  7
)


# ============================================================
# 7. MR-colocalization evidence landscape
# Strongly recommended for main/supplementary figure
# ============================================================

p_mr_coloc_landscape <- ev_plot %>%
  filter(
    mr_selected,
    coloc_success,
    is.finite(PP.H4),
    is.finite(pval)
  ) %>%
  ggplot(
    aes(
      x = log10P,
      y = PP.H4,
      color = Resource,
      size = conditional_pct_score
    )
  ) +
  geom_point(
    alpha = 0.65
  ) +
  geom_hline(
    yintercept = moderate_coloc_cutoff,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = strong_coloc_cutoff,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_size_continuous(
    range = c(1.5, 6)
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = expression(-log[10](P[MR])),
    y = "PP.H4",
    color = "Resource",
    size = "Conditional\npercentile score"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

save_plot_both(
  p_mr_coloc_landscape,
  file.path(dirs$qc, "09_MR_coloc_evidence_landscape"),
  8.5,
  6.5
)


# ============================================================
# 8. Evidence-class composition by resource
# Weak / Moderate / Strong
# ============================================================

evidence_composition <- ev_plot %>%
  filter(
    mr_selected,
    coloc_success
  ) %>%
  mutate(
    Evidence_class = case_when(
      PP.H4 >= strong_coloc_cutoff ~ "Strong",
      PP.H4 >= moderate_coloc_cutoff ~ "Moderate",
      TRUE ~ "Weak"
    ),
    Evidence_class = factor(
      Evidence_class,
      levels = c("Weak", "Moderate", "Strong")
    )
  ) %>%
  count(
    Resource,
    Evidence_class
  ) %>%
  group_by(Resource) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

p_evidence_composition <- ggplot(
  evidence_composition,
  aes(
    x = Resource,
    y = proportion,
    fill = Evidence_class
  )
) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(
      label = percent(proportion, accuracy = 0.1)
    ),
    position = position_stack(vjust = 0.5),
    size = 3.3
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Proportion of successfully colocalized pairs",
    fill = "Evidence"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_plot_both(
  p_evidence_composition,
  file.path(dirs$qc, "10_coloc_evidence_composition"),
  7,
  6
)


# ============================================================
# 9. Pipeline attrition as proportion of initial MR pairs
# More intuitive than raw counts alone
# ============================================================

pipeline_prop <- pipeline_qc %>%
  group_by(source) %>%
  mutate(
    baseline_n = n[stage == "MR_pairs"][1],
    proportion = n / baseline_n,
    Resource = recode(source, !!!source_labels)
  ) %>%
  ungroup()

p_pipeline_prop <- ggplot(
  pipeline_prop,
  aes(
    x = stage,
    y = proportion,
    group = Resource,
    color = Resource
  )
) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_text(
    aes(
      label = percent(proportion, accuracy = 0.1)
    ),
    vjust = -0.7,
    size = 3
  ) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, 1.08)
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Proportion of all MR-tested gene-context pairs",
    color = "Resource"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    panel.grid.minor = element_blank()
  )

save_plot_both(
  p_pipeline_prop,
  file.path(dirs$qc, "11_pipeline_attrition_proportion"),
  9,
  6
)


# ============================================================
# 10. Strong-support breadth across contexts
# Number of strong contexts per gene
# ============================================================

gene_strong_breadth <- ev_plot %>%
  group_by(
    gene_symbol,
    Resource
  ) %>%
  summarise(
    n_tested_contexts = n(),
    n_selected_contexts = sum(mr_selected, na.rm = TRUE),
    n_coloc_contexts = sum(coloc_success, na.rm = TRUE),
    n_strong_contexts = sum(
      strong_support,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_strong_contexts > 0
  )

p_strong_breadth <- gene_strong_breadth %>%
  ggplot(
    aes(
      x = n_strong_contexts,
      fill = Resource
    )
  ) +
  geom_histogram(
    binwidth = 1,
    boundary = 0.5,
    alpha = 0.6,
    position = "identity"
  ) +
  facet_wrap(
    ~Resource,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks()
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = paste0(
      "Number of contexts with strong support (PP.H4 ≥ ",
      strong_coloc_cutoff,
      ")"
    ),
    y = "Number of genes",
    fill = "Resource"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

save_plot_both(
  p_strong_breadth,
  file.path(dirs$qc, "12_gene_strong_context_breadth"),
  8,
  7
)


# ============================================================
# 11. Optional combined QC panel
# Requires patchwork, already loaded in main script
# ============================================================

p_qc_combined <-
  p_pipeline_prop /
  (p_raw_dist_color | p_evidence_composition) /
  p_mr_coloc_landscape +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "DynaMac-GP quality-control and evidence landscape"
  )

ggsave(
  file.path(
    dirs$qc,
    "QC_combined_publication_figure.pdf"
  ),
  p_qc_combined,
  width = 14,
  height = 16
)

ggsave(
  file.path(
    dirs$qc,
    "QC_combined_publication_figure.png"
  ),
  p_qc_combined,
  width = 14,
  height = 16,
  dpi = 600
)















