#!/usr/bin/env bash
# =============================================================================
# make_dummy_bed.sh
#
# Generates a PLINK binary (.bed/.bim/.fam) dummy dataset for testing the
# GMMAT gene-based analysis pipeline.
#
# Prerequisites: plink (1.9.x) must be in PATH
# Usage:  bash make_dummy_bed.sh
# =============================================================================

set -euo pipefail

OUTDIR="dummy_data"
PREFIX="${OUTDIR}/genotypes"

echo "Generating dummy .ped file..."

# Build a minimal PED from the FAM file + random genotypes.
# PED format: FID IID PAT MAT SEX PHENO geno1_A1 geno1_A2 geno2_A1 ...
# We generate random dosage 0/1/2 → homref / het / homalt calls.

Rscript - <<'RSCRIPT'
set.seed(99)
fam  <- read.table("dummy_data/genotypes.fam",  header=FALSE)
bim  <- read.table("dummy_data/genotypes.bim",  header=FALSE)
n    <- nrow(fam)
n_snp <- nrow(bim)

# Simulate genotypes with realistic MAF distribution
maf_vec <- c(runif(n_snp * 0.8, 0.05, 0.45),   # common
             runif(n_snp * 0.2, 0.001, 0.009))  # rare
maf_vec <- sample(maf_vec)[seq_len(n_snp)]

allele_pairs <- character(n_snp * n * 2)
idx <- 1L
for (j in seq_len(n_snp)) {
  maf  <- maf_vec[j]
  ref  <- as.character(bim$V5[j])
  alt  <- as.character(bim$V6[j])
  dos  <- sample(0:2, n, replace=TRUE,
                 prob=c((1-maf)^2, 2*maf*(1-maf), maf^2))
  a1 <- ifelse(dos == 0, ref, alt)
  a2 <- ifelse(dos == 2, alt, ref)
  allele_pairs[((j-1)*n*2 + 1):((j)*n*2)] <- c(rbind(a1, a2))
}

geno_mat <- matrix(allele_pairs, nrow=n, ncol=n_snp*2)
ped <- cbind(as.matrix(fam), geno_mat)

write.table(ped, "dummy_data/genotypes.ped",
            sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("PED written:", n, "samples x", n_snp, "variants\n")
RSCRIPT

echo "Converting PED → BED with PLINK..."
plink \
  --ped "${PREFIX}.ped" \
  --map "${PREFIX}.bim" \
  --make-bed \
  --out "${PREFIX}" \
  --allow-no-sex \
  --no-pheno

echo ""
echo "Done. PLINK binary files:"
ls -lh "${PREFIX}".{bed,bim,fam}
echo ""
echo "You can now run:"
echo "  Rscript 01_GMMAT_Common_Variant_Gene_Based.R"
echo "  Rscript 02_GMMAT_Rare_Variant_Gene_Based.R"
