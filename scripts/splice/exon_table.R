suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

x <- read_tsv(input, comment = "#", col_types = "c-cii-c-c",
              col_names = c("chr", "type", "start", "end", "strand", "info")) |>
    filter(type == "exon",
           chr %in% str_glue("chr{1:20}")) |>
    mutate(gene_id = str_match(info, 'gene_id "([^"]+)"')[, 2]) |>
    select(-type, -info)

write_tsv(x, output)
