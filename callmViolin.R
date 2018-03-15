#!/usr/bin/env Rscript

libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}
library(argparse)
options(scipen=999)


parser <- ArgumentParser()
parser$add_argument("--hichip", type = "character", required = FALSE,
  help = "A query loop file with columns being locus 1 and 2, genes 1 and 2 and optional columns.")
parser$add_argument("--din", type = "character", required = FALSE,
  help = "Feature file directory.")
parser$add_argument("--dout", type = "character", required = FALSE,
  help = "Pdffout directory.")
args <- parser$parse_args()
attach(args)

message("constructing loop and fet objects")
din <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/Spec_H3K27AC_ESC"
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = FALSE)
fet_fs <- dir(din, "*_heatmap.mat", recursive = TRUE, full.names = TRUE)
fet.obj <- fetConst(fet_fs, small = 0.05) 
mViolinPlot(loop.obj, fet.obj, dout = gsub("Spec_H3K27AC_ESC", "conFet", din))
