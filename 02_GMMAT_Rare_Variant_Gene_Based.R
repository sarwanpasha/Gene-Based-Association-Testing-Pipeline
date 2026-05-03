#!/usr/bin/env Rscript
# =============================================================================
# SCRIPT:  02_GMMAT_Rare_Variant_Gene_Based.R
# PURPOSE: Gene-based rare variant association testing using GMMAT SMMAT
#          (SKAT-O / burden)
# VERSION: 1.0.0
# =============================================================================
#
# DESCRIPTION:
#   Performs gene-based association testing for rare variants (MAF < 1%) using
#   GMMAT's SMMAT framework. Fits a null GLMM, builds a group file from a gene
#   annotation BED using Beta(1,25) SKAT weights, runs SMMAT with SKAT-O and
#   burden tests, and produces annotated results with Manhattan and QQ plots.
#
# USAGE:
#   Rscript 02_GMMAT_Rare_Variant_Gene_Based.R
#   -- or edit USER SETTINGS below and source interactively --
#
# REQUIREMENTS:
#   install.packages(c("data.table","ggplot2","ggrepel","dplyr","scales","doMC"))
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("GMMAT","SeqArray","SeqVarTools"))
#
# INPUT FILES:
#   1. PLINK binary set  : <PLINK_PREFIX>.{bed,bim,fam}
#   2. Phenotype file    : tab-separated, must contain IID, phenotype, covariates
#   3. Gene annotation   : tab-separated BED (chr, start, end, gene_name) – no header
#
# OUTPUT FILES (written to OUT_DIR):
#   smmat_rare_final.txt              – Raw SMMAT gene-level results
#   rare_variant_gene_results.txt     – Annotated results with gene coordinates
#   rare_variant_summary_stats.txt    – Per-method summary (lambda, n_sig, etc.)
#   rare_variant_manhattan_input.txt
#   rare_variant_manhattan.png
#   rare_variant_qqplot.png
#
# NOTES:
#   - SMMAT requires annotation/id (chr:pos:ref:alt) in the group file, NOT
#     integer variant.id from SeqArray — ensure your BIM/GDS is formatted so
#     the annotation/id field is populated (seqBED2GDS does this automatically).
#   - GDS must be closed before calling SMMAT() — it opens the file internally.
#   - When ncores > 1, per-thread score/var files are combined by SMMAT.meta().
#   - doMC is used for multi-threading (Linux/Mac). On Windows, set NTHREADS=1.
# =============================================================================

cat("\n========================================================\n")
cat("  GMMAT Gene-Based Rare Variant Analysis (SMMAT)\n")
cat("========================================================\n\n")

# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(GMMAT)
  library(SeqArray)
  library(SeqVarTools)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(scales)
  library(doMC)
})

# ── 1. USER SETTINGS ──────────────────────────────────────────────────────────
# ---- Paths -------------------------------------------------------------------
PLINK_PREFIX <- "/path/to/your/plink/genotype_file"   # no extension
PHENO_FILE   <- "/path/to/your/phenotype.txt"
GENE_BED     <- "/path/to/your/gene_annotation.bed"
OUT_DIR      <- "rare_variant_results"
GDS_FILE     <- file.path(OUT_DIR, "genotypes_rare.gds")

# If a GDS was already created by the common variant script, reuse it:
# COMMON_GDS   <- "common_variant_results/genotypes_common.gds"

# ---- Phenotype / model -------------------------------------------------------
PHENOTYPE  <- "OUTCOME"                          # column name of binary outcome (0/1)
COVARIATES <- c("Age", "Sex", "PC1", "PC2", "PC3")  # adjust as needed
FAMILY     <- binomial(link = "logit")           # binomial | gaussian

# ---- Variant filters ---------------------------------------------------------
MAF_RARE_MAX <- 0.01   # rare variant MAF upper bound
MAC_MIN      <- 3      # minimum allele count (match your QC threshold)

# ---- Weighting ---------------------------------------------------------------
# Beta(1,25): standard SKAT weights — heavily upweights ultra-rare variants
# Change to c(0.5,0.5) for flat weights (equal contribution per variant)
WEIGHTS_BETA <- c(1, 25)
MIN_RARE_VAR <- 2      # min rare variants required per gene to test

# ---- Computational -----------------------------------------------------------
NTHREADS    <- 4
RANDOM_SEED <- 42

# ---- Plot thresholds ---------------------------------------------------------
SIG_LINE     <- 2.5e-6   # gene-level Bonferroni (~0.05 / 20,000 genes)
SUGG_LINE    <- 1e-4     # suggestive threshold
TOP_N_LABELS <- 15

# ── 2. SETUP ──────────────────────────────────────────────────────────────────
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
registerDoMC(cores = NTHREADS)
set.seed(RANDOM_SEED)
cat("Output directory:", OUT_DIR, "\n\n")

# ── 3. LOAD & CLEAN PHENOTYPE DATA ────────────────────────────────────────────
cat("Loading phenotype data...\n")
pheno <- fread(PHENO_FILE, header = TRUE, sep = "\t",
               na.strings = c("NA", "", "."))
pheno <- as.data.frame(pheno)
cat("  Samples in phenotype file:", nrow(pheno), "\n")
cat("  Columns:", paste(names(pheno), collapse = ", "), "\n")

# Recode binary phenotype 1/2 → 0/1 if needed (PLINK convention)
if (PHENOTYPE %in% names(pheno)) {
  ph_vals <- unique(na.omit(pheno[[PHENOTYPE]]))
  if (all(ph_vals %in% c(1, 2))) {
    cat("  Recoding phenotype 1/2 → 0/1\n")
    pheno[[PHENOTYPE]] <- pheno[[PHENOTYPE]] - 1
  }
}

pheno$IID <- as.character(pheno$IID)
if (!"FID" %in% names(pheno)) pheno$FID <- pheno$IID
pheno$FID <- as.character(pheno$FID)

keep_cols   <- unique(c("FID", "IID", PHENOTYPE, COVARIATES))
keep_cols   <- keep_cols[keep_cols %in% names(pheno)]
pheno_clean <- pheno[, keep_cols, drop = FALSE]
pheno_clean <- pheno_clean[complete.cases(pheno_clean), ]
cat("  Samples after removing missing data:", nrow(pheno_clean), "\n")
cat("  Case/Control distribution:\n")
print(table(pheno_clean[[PHENOTYPE]]))
cat("\n")

# ── 4. PLINK → GDS CONVERSION ─────────────────────────────────────────────────
cat("Converting PLINK to GDS (or reusing existing)...\n")

# Reuse GDS from common variant script if available
common_gds <- "common_variant_results/genotypes_common.gds"
if (!file.exists(GDS_FILE) && file.exists(common_gds)) {
  cat("  Reusing GDS from common variant analysis:", common_gds, "\n")
  GDS_FILE <- common_gds
} else if (!file.exists(GDS_FILE)) {
  SeqArray::seqBED2GDS(
    bed.fn        = paste0(PLINK_PREFIX, ".bed"),
    fam.fn        = paste0(PLINK_PREFIX, ".fam"),
    bim.fn        = paste0(PLINK_PREFIX, ".bim"),
    out.gdsfn     = GDS_FILE,
    compress.geno = "ZIP_RA",
    verbose       = TRUE
  )
  cat("  GDS created:", GDS_FILE, "\n\n")
} else {
  cat("  GDS already exists – skipping conversion.\n\n")
}

# ── 5. FILTER RARE VARIANTS ───────────────────────────────────────────────────
cat("Filtering for rare variants (MAF <", MAF_RARE_MAX, ", MAC >=", MAC_MIN, ")...\n")
gds <- seqOpen(GDS_FILE)

gds_samples    <- seqGetData(gds, "sample.id")
common_samples <- intersect(gds_samples, pheno_clean$IID)
cat("  GDS samples:", length(gds_samples), "\n")
cat("  Overlapping with phenotype:", length(common_samples), "\n")

if (length(common_samples) < 10) {
  seqClose(gds)
  stop("Too few overlapping samples. Check IID matching between GDS and phenotype file.")
}

seqSetFilter(gds, sample.id = common_samples, verbose = FALSE)

af_all      <- seqAlleleFreq(gds, ref.allele = 0L)
maf_all     <- pmin(af_all, 1 - af_all)
variant_ids <- seqGetData(gds, "variant.id")
n_samples   <- length(common_samples)
mac_all     <- round(maf_all * 2 * n_samples)

rare_variants <- variant_ids[maf_all < MAF_RARE_MAX & mac_all >= MAC_MIN]
cat("  Total variants:", length(variant_ids), "\n")
cat("  Rare variants:", length(rare_variants), "\n\n")

seqSetFilter(gds, variant.id = rare_variants,
             sample.id = common_samples, verbose = FALSE)

maf_rare     <- pmin(seqAlleleFreq(gds, ref.allele = 0L), 1 - seqAlleleFreq(gds, ref.allele = 0L))
var_chr_rare <- seqGetData(gds, "chromosome")
var_pos_rare <- seqGetData(gds, "position")
var_ids_rare <- seqGetData(gds, "variant.id")

# annotation/id must be chr:pos:ref:alt — SMMAT matches variants via this field
var_annot_rare <- seqGetData(gds, "annotation/id")
cat("  Annotation ID sample (first 3):",
    paste(head(var_annot_rare, 3), collapse = ", "), "\n\n")

# ── 6. GENE ANNOTATION ────────────────────────────────────────────────────────
cat("Loading gene annotation...\n")
if (!file.exists(GENE_BED)) {
  seqClose(gds); stop("Gene annotation BED not found: ", GENE_BED)
}

gene_annot <- fread(GENE_BED, header = FALSE, sep = "\t")
if (ncol(gene_annot) < 4)
  stop("Gene BED must have ≥4 columns: chr, start, end, gene_name")
names(gene_annot)[1:4] <- c("chr", "start", "end", "gene")
gene_annot$chr <- sub("^chr", "", gene_annot$chr)
gene_annot     <- gene_annot[!duplicated(gene_annot$gene), ]
cat("  Genes loaded:", nrow(gene_annot), "\n\n")

# ── 7. MAP VARIANTS → GENES ───────────────────────────────────────────────────
cat("Mapping rare variants to genes...\n")

var_df_rare <- data.frame(
  variant.id = var_ids_rare,
  annot.id   = var_annot_rare,
  chr        = as.character(var_chr_rare),
  pos        = var_pos_rare,
  maf        = maf_rare,
  stringsAsFactors = FALSE
)

gene_list_rare <- list()
gene_maf_list  <- list()

for (i in seq_len(nrow(gene_annot))) {
  g_chr   <- as.character(gene_annot$chr[i])
  g_start <- gene_annot$start[i]
  g_end   <- gene_annot$end[i]
  g_name  <- gene_annot$gene[i]
  idx     <- which(var_df_rare$chr == g_chr &
                   var_df_rare$pos >= g_start &
                   var_df_rare$pos <= g_end)
  if (length(idx) >= MIN_RARE_VAR) {
    gene_list_rare[[g_name]] <- var_df_rare$annot.id[idx]
    gene_maf_list[[g_name]]  <- var_df_rare$maf[idx]
  }
}

rare_counts <- sapply(gene_list_rare, length)
cat("  Genes with >=", MIN_RARE_VAR, "rare variants:", length(gene_list_rare), "\n")
cat("  Variants/gene — median:", median(rare_counts),
    "  range: [", min(rare_counts), "-", max(rare_counts), "]\n\n")

if (length(gene_list_rare) == 0) {
  seqClose(gds)
  stop("No genes found with rare variants. Check chromosome naming and BED coordinates.")
}

# ── 8. NULL MODEL ─────────────────────────────────────────────────────────────
cat("Fitting null GLMM...\n")

pheno_aligned <- pheno_clean[match(common_samples, pheno_clean$IID), ]
rownames(pheno_aligned) <- pheno_aligned$IID

valid_covars <- COVARIATES[COVARIATES %in% names(pheno_aligned)]
covar_str    <- paste(valid_covars, collapse = " + ")
formula_str  <- if (nchar(covar_str) > 0) {
  paste(PHENOTYPE, "~", covar_str)
} else {
  paste(PHENOTYPE, "~ 1")
}
cat("  Formula:", formula_str, "\n")

null_model <- GMMAT::glmmkin(
  fixed        = as.formula(formula_str),
  data         = pheno_aligned,
  kins         = NULL,
  id           = "IID",
  family       = FAMILY,
  method       = "REML",
  method.optim = "AI"
)
cat("  Null model fitted successfully.\n\n")

# ── 9. BUILD GROUP FILE ───────────────────────────────────────────────────────
# Format (no header): group | chr | pos | ref | alt | weight
# Beta(1,25) weights strongly upweight ultra-rare variants (standard for SKAT).
cat("Building SMMAT group file...\n")

compute_beta_weights <- function(maf_vec, a1 = WEIGHTS_BETA[1], a2 = WEIGHTS_BETA[2]) {
  w  <- dbeta(maf_vec, shape1 = a1, shape2 = a2)
  mx <- max(w, na.rm = TRUE)
  if (mx == 0) return(rep(1, length(w)))
  w / mx
}

group_rows <- vector("list", sum(sapply(gene_list_rare, length)))
idx_fill   <- 1L

for (gene_name in names(gene_list_rare)) {
  annot_ids <- gene_list_rare[[gene_name]]
  gene_mafs <- gene_maf_list[[gene_name]]
  beta_wts  <- compute_beta_weights(gene_mafs)
  for (j in seq_along(annot_ids)) {
    parts <- strsplit(annot_ids[j], ":")[[1]]
    group_rows[[idx_fill]] <- data.frame(
      group  = gene_name,
      chr    = parts[1], pos = parts[2],
      ref    = parts[3], alt = parts[4],
      weight = round(beta_wts[j], 6),
      stringsAsFactors = FALSE
    )
    idx_fill <- idx_fill + 1L
  }
}

group_df   <- do.call(rbind, group_rows)
group_file <- file.path(OUT_DIR, "gene_groups_rare.txt")
fwrite(group_df, group_file, sep = "\t", col.names = FALSE, quote = FALSE, eol = "\n")
cat("  Group file:", nrow(group_df), "rows,", length(gene_list_rare), "genes\n")
cat("  First 3 rows:\n"); print(head(group_df, 3)); cat("\n")

# ── 10. RUN SMMAT ─────────────────────────────────────────────────────────────
seqClose(gds)
cat("GDS closed. Running SMMAT...\n")
cat("  Tests: Efficient SKAT-O (E), Burden (B)\n")
cat("  Genes:", length(gene_list_rare), "\n\n")

GMMAT::SMMAT(
  null.obj         = null_model,
  geno.file        = GDS_FILE,
  group.file       = group_file,
  group.file.sep   = "\t",
  meta.file.prefix = file.path(OUT_DIR, "smmat_rare"),
  MAF.range        = c(1e-7, MAF_RARE_MAX),
  MAF.weights.beta = WEIGHTS_BETA,
  miss.cutoff      = 0.15,
  missing.method   = "impute2mean",
  method           = "davies",
  tests            = c("E", "B"),
  rho              = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
  use.minor.allele = TRUE,
  auto.flip        = TRUE,
  ncores           = NTHREADS,
  verbose          = TRUE
)

cat("Combining per-thread results via SMMAT.meta()...\n")
smmat_res <- GMMAT::SMMAT.meta(
  meta.files.prefix = file.path(OUT_DIR, "smmat_rare"),
  n.files           = NTHREADS,
  group.file        = group_file,
  group.file.sep    = "\t",
  MAF.range         = c(1e-7, MAF_RARE_MAX),
  MAF.weights.beta  = WEIGHTS_BETA,
  miss.cutoff       = 0.15,
  method            = "davies",
  tests             = c("E", "B"),
  rho               = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
  use.minor.allele  = TRUE,
  verbose           = TRUE
)

out_smmat <- file.path(OUT_DIR, "smmat_rare_final.txt")
fwrite(smmat_res, out_smmat, sep = "\t", quote = FALSE)
cat("  SMMAT results written:", out_smmat, "\n\n")

# ── 11. ANNOTATE RESULTS ──────────────────────────────────────────────────────
smmat_res <- as.data.table(smmat_res)
if ("group" %in% names(smmat_res)) setnames(smmat_res, "group", "gene")

cat("  Genes with results:", nrow(smmat_res), "\n")
cat("  Columns:", paste(names(smmat_res), collapse = ", "), "\n\n")

gene_coords   <- gene_annot[, .(chr, start, end, gene)]
combined_rare <- merge(smmat_res, gene_coords, by = "gene", all.x = TRUE)
combined_rare$chr_num <- suppressWarnings(as.numeric(combined_rare$chr))
combined_rare <- combined_rare[!is.na(combined_rare$chr_num), ]
combined_rare <- combined_rare[order(combined_rare$chr_num, combined_rare$start), ]
if ("n.variants" %in% names(combined_rare))
  setnames(combined_rare, "n.variants", "n_var")

out_results <- file.path(OUT_DIR, "rare_variant_gene_results.txt")
fwrite(combined_rare, out_results, sep = "\t", quote = FALSE, na = "NA")
cat("Annotated results saved:", out_results, "\n\n")

# ── 12. SUMMARY STATISTICS ────────────────────────────────────────────────────
lambda_gc_fn <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 2) return(NA_real_)
  chi2 <- qchisq(pvals, df = 1, lower.tail = FALSE)
  round(median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1), 4)
}

p_cols <- intersect(
  c("p.E","p.B","p.S","p.O","E.pval","B.pval","S.pval","O.pval"),
  names(combined_rare)
)
cat("P-value columns found:", paste(p_cols, collapse = ", "), "\n")

summary_list <- lapply(p_cols, function(pc) {
  pv <- combined_rare[[pc]]
  n  <- sum(!is.na(pv))
  data.frame(
    method           = pc,
    n_genes_tested   = n,
    n_sig_bonferroni = sum(pv < (0.05 / n), na.rm = TRUE),
    n_sig_threshold  = sum(pv < SIG_LINE,   na.rm = TRUE),
    n_suggestive     = sum(pv < SUGG_LINE,  na.rm = TRUE),
    min_p            = min(pv,    na.rm = TRUE),
    median_p         = median(pv, na.rm = TRUE),
    lambda_gc        = lambda_gc_fn(pv),
    stringsAsFactors = FALSE
  )
})
summary_stats <- do.call(rbind, summary_list)

out_summary <- file.path(OUT_DIR, "rare_variant_summary_stats.txt")
fwrite(summary_stats, out_summary, sep = "\t", quote = FALSE)
cat("Summary statistics:\n"); print(as.data.frame(summary_stats)); cat("\n")

# ── 13. MANHATTAN PLOT ────────────────────────────────────────────────────────
primary_pcol <- if ("p.E"    %in% names(combined_rare)) "p.E"    else
                if ("E.pval" %in% names(combined_rare)) "E.pval" else p_cols[1]
cat("Primary p-value column for plots:", primary_pcol, "\n")

manhattan_data_rare <- combined_rare %>%
  as.data.frame() %>%
  mutate(
    CHR      = chr_num,
    BP       = as.integer((start + end) / 2),
    P        = .data[[primary_pcol]],
    SNP      = gene,
    n_rare   = if ("n_var"     %in% names(.)) n_var     else NA_integer_,
    mean_maf = if ("freq.mean" %in% names(.)) freq.mean else NA_real_
  ) %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P), P > 0, P <= 1) %>%
  mutate(negLogP = -log10(P)) %>%
  arrange(CHR, BP)

out_manhattan_input <- file.path(OUT_DIR, "rare_variant_manhattan_input.txt")
fwrite(
  manhattan_data_rare[, intersect(c("CHR","BP","SNP","P","negLogP","n_rare","mean_maf"),
                                  names(manhattan_data_rare))],
  out_manhattan_input, sep = "\t", quote = FALSE
)

lambda_gc_rare <- lambda_gc_fn(manhattan_data_rare$P)
cat("Lambda GC:", lambda_gc_rare, "\n\n")

cat("Generating Manhattan plot...\n")

chr_sizes_rare <- manhattan_data_rare %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(bp_add = cumsum(as.numeric(lag(max_bp, default = 0))) + (CHR - 1) * 5e7)

manhattan_data_rare <- manhattan_data_rare %>%
  left_join(chr_sizes_rare[, c("CHR","bp_add")], by = "CHR") %>%
  mutate(BP_cum = BP + bp_add)

chr_labels_rare <- manhattan_data_rare %>%
  group_by(CHR) %>% summarise(center = mean(BP_cum), .groups = "drop")

chr_colors_rare <- rep(c("#6A0572","#C77DFF"), 12)[seq_along(unique(manhattan_data_rare$CHR))]
names(chr_colors_rare) <- as.character(sort(unique(manhattan_data_rare$CHR)))

manhattan_data_rare <- manhattan_data_rare %>%
  mutate(
    color_group = as.character(CHR),
    point_size  = ifelse(P < SIG_LINE, 2.5, ifelse(P < SUGG_LINE, 1.8, 1.0)),
    is_sig      = P < SIG_LINE,
    is_sugg     = P < SUGG_LINE & P >= SIG_LINE
  )

top_genes_rare <- manhattan_data_rare %>%
  filter(P < SUGG_LINE) %>% arrange(P) %>% slice_head(n = TOP_N_LABELS)

p_manhattan_rare <- ggplot(manhattan_data_rare,
                           aes(x = BP_cum, y = negLogP, color = color_group)) +
  geom_point(aes(size = point_size), alpha = 0.7, shape = 16) +
  geom_point(data = filter(manhattan_data_rare, is_sugg),
             color = "#FF9F1C", size = 2.2, shape = 16) +
  geom_point(data = filter(manhattan_data_rare, is_sig),
             color = "#D62828", size = 3.2, shape = 17) +
  geom_hline(yintercept = -log10(SIG_LINE),  linetype = "dashed",
             color = "#D62828", linewidth = 0.9) +
  geom_hline(yintercept = -log10(SUGG_LINE), linetype = "dotted",
             color = "#FF9F1C", linewidth = 0.7) +
  ggrepel::geom_label_repel(
    data = top_genes_rare, aes(label = SNP),
    size = 2.8, fontface = "bold.italic", color = "black",
    fill = scales::alpha("white", 0.85),
    box.padding = 0.4, point.padding = 0.3,
    segment.color = "gray40", segment.size = 0.3,
    max.overlaps = 20, show.legend = FALSE
  ) +
  scale_color_manual(values = chr_colors_rare, guide = "none") +
  scale_size_identity() +
  scale_x_continuous(breaks = chr_labels_rare$center, labels = chr_labels_rare$CHR,
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(
    title    = paste0("Rare Variant Gene-Based Analysis (SMMAT-",
                      sub("p\\.", "", primary_pcol), ")"),
    subtitle = bquote(lambda[GC] == .(lambda_gc_rare) ~
                        " | MAF <" ~ .(MAF_RARE_MAX) ~ ", MAC \u2265" ~ .(MAC_MIN) ~
                        " | n =" ~ .(nrow(pheno_aligned)) ~ "samples"),
    x = "Chromosome", y = expression(-log[10](p))
  ) +
  annotate("text", x = max(manhattan_data_rare$BP_cum) * 0.02,
           y = -log10(SIG_LINE)  + 0.15,
           label = paste0("p = ", format(SIG_LINE,  scientific = TRUE)),
           color = "#D62828", size = 2.8, hjust = 0) +
  annotate("text", x = max(manhattan_data_rare$BP_cum) * 0.02,
           y = -log10(SUGG_LINE) + 0.15,
           label = paste0("p = ", format(SUGG_LINE, scientific = TRUE)),
           color = "#FF9F1C", size = 2.8, hjust = 0) +
  theme_bw(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray30"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "gray40"),
    axis.title         = element_text(face = "bold"),
    axis.text.x        = element_text(size = 9),
    plot.margin        = margin(10, 20, 10, 10)
  )

out_manhattan_rare <- file.path(OUT_DIR, "rare_variant_manhattan.png")
ggsave(out_manhattan_rare, p_manhattan_rare, width = 16, height = 6.5,
       dpi = 300, bg = "white")
cat("  Manhattan plot saved:", out_manhattan_rare, "\n")

# ── 14. QQ PLOT ───────────────────────────────────────────────────────────────
cat("Generating QQ plot...\n")

qq_data_rare <- manhattan_data_rare %>%
  filter(!is.na(P), P > 0, P <= 1) %>%
  arrange(P) %>%
  mutate(
    observed  = -log10(P),
    expected  = -log10(ppoints(n())),
    color_cat = case_when(
      P < SIG_LINE  ~ "Significant",
      P < SUGG_LINE ~ "Suggestive",
      TRUE          ~ "Nominal"
    )
  )

n_tests_rare <- nrow(qq_data_rare)
qq_data_rare <- qq_data_rare %>%
  mutate(
    ci_upper = sort(-log10(qbeta(0.025, seq_len(n_tests_rare), n_tests_rare:1)),
                   decreasing = TRUE),
    ci_lower = sort(-log10(qbeta(0.975, seq_len(n_tests_rare), n_tests_rare:1)),
                   decreasing = TRUE)
  )

p_qq_rare <- ggplot(qq_data_rare, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "gray85", alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "gray40", linewidth = 0.8,
              linetype = "dashed") +
  geom_point(aes(color = color_cat), size = 1.8, alpha = 0.8, shape = 16) +
  geom_point(data = filter(qq_data_rare, color_cat != "Nominal"),
             aes(color = color_cat), size = 2.8, shape = 16) +
  ggrepel::geom_label_repel(
    data = filter(qq_data_rare, P < SUGG_LINE) %>% slice_head(n = TOP_N_LABELS),
    aes(label = SNP, color = color_cat),
    size = 2.5, fontface = "bold.italic",
    fill = scales::alpha("white", 0.85),
    box.padding = 0.4, show.legend = FALSE, max.overlaps = 15
  ) +
  scale_color_manual(
    values = c("Significant" = "#D62828", "Suggestive" = "#FF9F1C",
               "Nominal" = "#6A0572"),
    name = "Significance", breaks = c("Significant", "Suggestive", "Nominal")
  ) +
  labs(
    title    = paste0("QQ Plot – Rare Variant Gene-Based Analysis (SMMAT-",
                      sub("p\\.", "", primary_pcol), ")"),
    subtitle = bquote(lambda[GC] == .(lambda_gc_rare) ~
                        " | MAF <" ~ .(MAF_RARE_MAX)),
    x = expression("Expected" ~ -log[10](p)),
    y = expression("Observed"  ~ -log[10](p))
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "gray30"),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "gray40"),
    legend.position   = c(0.15, 0.85),
    legend.background = element_rect(fill = scales::alpha("white", 0.8),
                                     color = "gray70"),
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9),
    axis.title        = element_text(face = "bold")
  ) +
  annotate("label", x = Inf, y = -Inf,
           label = paste0("\u03bb = ", lambda_gc_rare),
           hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold",
           fill = "#F3E5FF", color = "#4A0066", label.size = 0.6)

out_qq_rare <- file.path(OUT_DIR, "rare_variant_qqplot.png")
ggsave(out_qq_rare, p_qq_rare, width = 8, height = 8, dpi = 300, bg = "white")
cat("  QQ plot saved:", out_qq_rare, "\n")

# ── 15. TOP RESULTS TABLE ─────────────────────────────────────────────────────
cat("\nTop 20 genes by", primary_pcol, ":\n")
top20 <- combined_rare %>%
  as.data.frame() %>%
  arrange(.data[[primary_pcol]]) %>%
  select(gene, chr, any_of(c("n_var","freq.mean")), any_of(p_cols)) %>%
  head(20)
print(top20, digits = 4)

# ── 16. DONE ──────────────────────────────────────────────────────────────────
# Note: GDS was closed before SMMAT — do NOT call seqClose(gds) again
cat("\n========================================================\n")
cat("  RARE VARIANT ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("  Output files:\n")
cat("   \u2022 SMMAT results      :", out_smmat, "\n")
cat("   \u2022 Annotated results  :", out_results, "\n")
cat("   \u2022 Summary statistics :", out_summary, "\n")
cat("   \u2022 Manhattan input    :", out_manhattan_input, "\n")
cat("   \u2022 Manhattan plot     :", out_manhattan_rare, "\n")
cat("   \u2022 QQ plot            :", out_qq_rare, "\n")
cat("  Lambda GC (", primary_pcol, "):", lambda_gc_rare, "\n\n")
cat("Session info:\n"); print(sessionInfo())
