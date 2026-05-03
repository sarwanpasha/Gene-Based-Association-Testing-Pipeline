# GMMAT Gene-Based Association Testing Pipeline

A two-script R pipeline for gene-based association testing of **common** and **rare** variants using [GMMAT](https://github.com/hanchenphd/GMMAT) SMMAT (Sequence Kernel Association Meta-analysis Test). Designed for binary outcomes (e.g., disease case/control) from PLINK-formatted whole-genome or whole-exome sequencing data.

---

## Contents

```
.
├── 01_GMMAT_Common_Variant_Gene_Based.R   # Common variant analysis (MAF ≥ 1%)
├── 02_GMMAT_Rare_Variant_Gene_Based.R     # Rare variant analysis   (MAF < 1%)
└── dummy_data/
    ├── phenotype.txt                      # 200 simulated samples
    ├── genotypes.bim                      # 300 SNPs across 3 chromosomes
    ├── genotypes.fam                      # Sample metadata
    ├── gene_annotation.bed                # 10 genes across 3 chromosomes
    ├── genotype_preview.csv               # Preview matrix (not used by scripts)
    └── make_dummy_bed.sh                  # Generates .bed file from .bim/.fam
```

---

## What Each Script Does

### `01_GMMAT_Common_Variant_Gene_Based.R`

1. Loads and cleans the phenotype/covariate file
2. Converts PLINK binary to GDS format (SeqArray)
3. Filters variants to MAF ≥ `MAF_COMMON` (default 1%) and MAC ≥ `MAC_MIN`
4. Maps filtered variants to gene windows from a BED annotation file
5. Fits a null GLMM (`glmmkin`) with user-specified covariates
6. Builds a SMMAT group file using flat Beta(0.5, 0.5) weights (uniform weighting)
7. Runs `SMMAT()` with SKAT (`S`), burden (`B`), and SKAT-O (`O`) tests
8. Combines per-thread results with `SMMAT.meta()`
9. Writes annotated results, summary statistics, Manhattan plot, and QQ plot

### `02_GMMAT_Rare_Variant_Gene_Based.R`

Same workflow as above, but:
- Filters to MAF < `MAF_RARE_MAX` (default 1%) with MAC ≥ `MAC_MIN`
- Uses Beta(1, 25) weights (SKAT standard — strongly upweights ultra-rare variants)
- Runs efficient SKAT-O (`E`) and burden (`B`) tests
- Reuses the GDS file from script 01 if available (avoids re-conversion)

---

## Requirements

### R packages

```r
# CRAN
install.packages(c("data.table", "ggplot2", "ggrepel", "dplyr", "scales", "doMC"))

# Bioconductor
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GMMAT", "SeqArray", "SeqVarTools"))
```

**Tested with:** R ≥ 4.2, GMMAT ≥ 1.4.1, SeqArray ≥ 1.38.0

> **Windows users:** `doMC` is not available on Windows. Set `NTHREADS <- 1` in the USER SETTINGS block of each script.

### External tools

- [PLINK 1.9](https://www.cog-genomics.org/plink/) — required only to generate the dummy BED file via `make_dummy_bed.sh`; not needed to run the R scripts themselves

---

## Input File Formats

### 1. PLINK binary set (`.bed` / `.bim` / `.fam`)

Standard PLINK 1.9 binary format. The BIM SNP ID column should be in `chr:pos:ref:alt` format — this is used as the `annotation/id` field in GDS and is how SMMAT matches variants in the group file.

If your BIM file uses `rs` IDs, update the BIM `$V2` column before conversion:
```r
bim <- read.table("genotypes.bim", header=FALSE)
bim$V2 <- paste(bim$V1, bim$V4, bim$V5, bim$V6, sep=":")
write.table(bim, "genotypes.bim", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

### 2. Phenotype file (tab-separated, with header)

Must contain at minimum: `IID`, your phenotype column, and all covariate columns.

| FID | IID | OUTCOME | Age | Sex | PC1 | PC2 | PC3 |
|---|---|---|---|---|---|---|---|
| FAM001 | SAMP001 | 1 | 72.3 | 0 | 0.123 | -0.045 | 0.891 |
| FAM002 | SAMP002 | 0 | 65.1 | 1 | -0.234 | 0.312 | -0.456 |

- Binary phenotype should be coded `0` (control) / `1` (case). If coded `1`/`2` (PLINK convention), the script auto-recodes.
- `FID` is optional; if absent, `IID` is used for both columns.

### 3. Gene annotation BED (tab-separated, **no header**)

| col 1 | col 2 | col 3 | col 4 |
|---|---|---|---|
| chromosome | start (bp) | end (bp) | gene name |

```
1	1000	5000	GENE_A
1	50000	55000	GENE_B
2	5000	12000	GENE_C
```

- `chr` prefix is stripped automatically (`chr1` → `1`).
- Duplicate gene names are removed (first occurrence kept).
- Coordinates should match the genome build of your genotype data.

---

## Usage

### Step 1 – Generate dummy PLINK BED (first-time test only)

```bash
cd dummy_data
bash make_dummy_bed.sh
cd ..
```

This produces `dummy_data/genotypes.bed` using PLINK.

### Step 2 – Configure paths in each script

Open `01_GMMAT_Common_Variant_Gene_Based.R` and edit the **USER SETTINGS** block:

```r
PLINK_PREFIX <- "dummy_data/genotypes"         # path to your .bed/.bim/.fam
PHENO_FILE   <- "dummy_data/phenotype.txt"
GENE_BED     <- "dummy_data/gene_annotation.bed"
OUT_DIR      <- "common_variant_results"

PHENOTYPE    <- "OUTCOME"
COVARIATES   <- c("Age", "Sex", "PC1", "PC2", "PC3")
```

Repeat for `02_GMMAT_Rare_Variant_Gene_Based.R`.

### Step 3 – Run

```bash
Rscript 01_GMMAT_Common_Variant_Gene_Based.R
Rscript 02_GMMAT_Rare_Variant_Gene_Based.R
```

Or run sequentially in one step:

```bash
Rscript 01_GMMAT_Common_Variant_Gene_Based.R && \
Rscript 02_GMMAT_Rare_Variant_Gene_Based.R
```

---

## Output Files

Both scripts write to their respective `OUT_DIR` (default: `common_variant_results/` and `rare_variant_results/`).

| File | Description |
|---|---|
| `smmat_*_final.txt` | Raw SMMAT output: one row per gene, p-values for each test |
| `*_gene_results.txt` | Annotated results merged with gene coordinates |
| `*_summary_stats.txt` | Per-test summary: n genes, n significant, lambda GC, min p |
| `*_manhattan_input.txt` | Cleaned data ready for custom plotting |
| `*_manhattan.png` | Manhattan plot (gene-level, colored by chromosome) |
| `*_qqplot.png` | QQ plot with 95% CI ribbon and genomic inflation factor |

### SMMAT output columns

| Column | Description |
|---|---|
| `gene` | Gene name (from group file) |
| `n_var` | Number of variants in the gene passing filters |
| `freq.min/mean/max` | Min/mean/max allele frequency across variants in gene |
| `p.S` / `S.pval` | SKAT p-value |
| `p.B` / `B.pval` | Burden test p-value |
| `p.O` / `O.pval` | SKAT-O p-value (common) |
| `p.E` / `E.pval` | Efficient SKAT-O p-value (rare) |

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `MAF_COMMON` | `0.01` | Lower MAF bound for common variant analysis |
| `MAF_RARE_MAX` | `0.01` | Upper MAF bound for rare variant analysis |
| `MAC_MIN` | `3` | Minimum minor allele count (should match your QC threshold) |
| `WEIGHTS_BETA` | `c(0.5,0.5)` (common) / `c(1,25)` (rare) | Beta distribution parameters for variant weighting |
| `MIN_COMMON_VAR` | `2` | Minimum variants per gene to include in testing |
| `SIG_LINE` | `2.5e-6` | Significance threshold (~Bonferroni for 20,000 genes) |
| `SUGG_LINE` | `1e-4` | Suggestive significance threshold |
| `NTHREADS` | `4` | Number of parallel threads for SMMAT |

---

## Design Notes

**Weighting:** Common variant analysis uses flat Beta(0.5, 0.5) weights, giving equal contribution to each variant. Rare variant analysis uses Beta(1, 25), which strongly upweights ultra-rare variants — this is the standard SKAT weighting scheme. Change `WEIGHTS_BETA` to `c(1, 25)` in the common variant script if you prefer SKAT-style weights.

**GDS reuse:** Script 02 checks for the GDS generated by script 01 and reuses it if present, avoiding redundant conversion. If you use different PLINK files for each script, set `GDS_FILE` to separate paths.

**SKAT-O for rare variants:** The efficient SKAT-O (`test = "E"`) is used in script 02. This is more powerful than standard SKAT-O for rare variants with small sample sizes.

**Null model:** No kinship matrix (`kins = NULL`) is used by default. For related samples, pass a precomputed GRM or kinship matrix to `glmmkin()` via the `kins` argument.

---

## Citation

If you use this pipeline, please cite:

- **GMMAT / SMMAT:** Chen H, et al. (2019). *Efficient variant set mixed model association tests for continuous and binary traits in large-scale whole-genome sequencing studies.* Am J Hum Genet. 104(2):260–274. https://doi.org/10.1016/j.ajhg.2018.12.012

- **SeqArray:** Zheng X, et al. (2017). *SeqArray — A storage-efficient high-performance data format for WGS variant calls.* Bioinformatics. 33(15):2251–2257.

---

## License

MIT
