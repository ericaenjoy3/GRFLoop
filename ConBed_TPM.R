#!/usr/bin/env Rscript

###
# Bigwig signals under peaks that overlapped with large windows
###

libfs <- c("data.table", "tidyverse", "gtools", "argparse", "RNA")
invisible(sapply(libfs, function(f)
	suppressPackageStartupMessages(require(f, character.only = TRUE))))
options(scipen = 999)

parser <- ArgumentParser()
parser$add_argument("--windowf", type = "character", required = FALSE,
  help = "A bed file of large windows (> 1 kb) with the 6th column of gene id.")
parser$add_argument("--proteincoding",
  help = "Switch to use protein-coding genes only.")
parser$add_argument("--matf", type = "character", required = FALSE,
  help = "Output signal matrix file.")
args <- parser$parse_args()
attach(args)

tpmf <- path.expand("~/athena/RNA/RNA_seq/DF5154_2017_08_25/salmon/TPM_rd_merge.txt")

# link gene id with gene expression
dat <- fread(windowf, header = FALSE)
setnames(dat, c("chr", "start", "end", "loc", "type", "gene"))

tpm <- SepTPMCnt(tpmf)$tpm.grp
tpm <- data.table(gene = rownames(tpm), tpm, key = "gene")
gid <- copy(dat[["gene"]])

dat[, `:=`(ESC = tpm[gid, ESC])]

if (!is.null(proteincoding)) {
	dat[!grepl(".+\\|protein_coding$", gene), `:=`(ESC = as.numeric(NA))]
}

# mean gene expression with the window (chr, start, end), considering unique genes
mean_tpm <- copy(dat[, {idx = !duplicated(gene); mean(ESC[idx], na.rm = TRUE)}, by = .(chr, start, end)])
setkeyv(mean_tpm, c("chr", "start", "end"))

# repeat output of gene expression and duplicate gene expression for the same window
dat[, `:=`(TPM_ESC = copy(mean_tpm[dat[, .(chr, start, end)], V1])) ]

write.table(dat[, .(TPM_ESC)], file = matf, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


