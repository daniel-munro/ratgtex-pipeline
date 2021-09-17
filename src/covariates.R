args <- commandArgs(trailingOnly = TRUE)
GENO_SIM_FILE <- args[1]
BED_FILE <- args[2]
N_PCS <- as.integer(args[3])
OUT_FILE <- args[4]

df <- read.delim(BED_FILE, check.names = FALSE, row.names = "gene_id")[, -(1:3)]
df <- df[apply(df, 1, var) != 0, ]
pca <- prcomp(t(df), center = TRUE, scale = TRUE)
pcs <- round(pca$x[, 1:N_PCS], 6)

# pcs_t <- t(pcs)

pcs_t <- data.frame(ID = colnames(pcs), t(pcs))
colnames(pcs_t) <- c("ID", rownames(pcs))  # column names starting with digits get 'fixed'

geno <- read.delim(GENO_SIM_FILE, check.names = FALSE)
geno$ID <- paste("simil", sub("-N", "", geno$ID), sep = "_")
geno <- geno[, colnames(pcs_t)]
df <- rbind(geno, pcs_t)

write.table(df, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
