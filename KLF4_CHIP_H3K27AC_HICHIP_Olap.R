#!/usr/bin/env Rscript

###
# Classify H3K27AC Hi-ChIP anchors by KLF4 ChIP-seq peaks
# Filters:
# (1) Filter hichip chr by chr in chip-seq peaks
###

libfs <- c("data.table", "tidyverse", "gtools", "igraph", "RJSONIO", "ggplot2", "argparse")
invisible(sapply(libfs, function(f)suppressPackageStartupMessages(require(f, character.only = T))))
options(scipen = 999)

parser <- ArgumentParser()
parser$add_argument("--hichip", type = "character", required = FALSE,
  help = "A query loop file with columns being locus 1 and 2, genes 1 and 2 and optional columns.")
parser$add_argument("--chip", type = "character", nargs='+', required = FALSE,
  help = "Chip-seq peaks of KLF4.")
parser$add_argument("--bed", type = "character", nargs='+', required = FALSE,
  help = "Output bed file with the 4th column indicating overlap status.")
args <- parser$parse_args()
attach(args)

dat <- fread(hichip, header = TRUE, sep = "\t", na.strings=c("N/A", ""))
dat <- rbindlist(list(dat[, .(loc1Chr, loc1Start, loc1End)], dat[, .(loc2Chr, loc2Start, loc2End)]), use.names = FALSE) %>% unique()

# Chip-seq peaks
chip_dat <- fread(chip, header = FALSE)[, 1:3, with = FALSE]
setnames(chip_dat, c("chr", "start", "end"))
setkeyv(chip_dat, c("chr", "start", "end"))
chip_dat[, rowid := 1:nrow(chip_dat)]

# (1) Filter hichip chr by chr in chip-seq peaks
message(dat[!loc1Chr %in% chip_dat[, unique(chr)], which = TRUE], " loops removed due to chr not present in the chip-seq file.")
dat <- dat[loc1Chr %in% chip_dat[, unique(chr)],]

setkeyv(dat, c("loc1Chr", "loc1Start", "loc1End"))
idx <- foverlaps(dat, chip_dat, which = TRUE, nomatch = 0)$xid
dat[idx, status := "With_KLF4"]
dat[is.na(status), status := "No_KLF4"]

write.table(dat, file = bed, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

