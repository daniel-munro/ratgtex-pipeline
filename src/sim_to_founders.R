suppressPackageStartupMessages(library(VariantAnnotation))

load_geno <- function(filename) {
    gt <- readGT(filename)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2,
                                       "0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)[x])
    # rownames(geno) <- rownames(gt)
    # If ID is missing, ID includes alleles, e.g. 'chr1:4315_G/GT'
    rownames(geno) <- sapply(strsplit(rownames(gt), "_", fixed = TRUE), function(x) x[1])
    geno
}

args <- commandArgs(trailingOnly = TRUE)
POP_VCF <- args[1]
FOUNDER_VCF <- args[2]
OUTPUT <- args[3]

cat("Loading population genotypes...\n")
pop <- load_geno(POP_VCF)
cat("Loading founder genotypes...\n")
founder <- load_geno(FOUNDER_VCF)

pop_ref <- setNames(as.character(ref(readVcf(POP_VCF))), rownames(pop))
founder_ref <- setNames(as.character(ref(readVcf(FOUNDER_VCF))), rownames(founder))
snps <- intersect(rownames(pop), rownames(founder))
snps2 <- snps[pop_ref[snps] == founder_ref[snps]]
cat("Removing", length(snps) - length(snps2), "of", length(snps), "SNPs due to ref mismatch\n")
pop <- pop[snps2, ]
founder <- founder[snps2, ]
# pop_ref <- pop_ref[snps, ]
# founder_ref <- founder_ref[snps, ]

cat("Calculating similarities...\n")
sim <- matrix(NA, nrow = ncol(founder), ncol = ncol(pop), dimnames = list(colnames(founder), colnames(pop)))
for (i in 1:ncol(founder)) {
    for (j in 1:ncol(pop)) {
        sim[i, j] <- mean(2 - abs(founder[, i] - pop[, j]), na.rm = TRUE) / 2
    }
}
sim <- as.data.frame(sim)
df <- data.frame(ID = rownames(sim), sim)
colnames(df) <- c("ID", colnames(sim))  # column names starting with digits get 'fixed'

write.table(df, OUTPUT, sep = "\t", quote = FALSE, row.names = FALSE)
