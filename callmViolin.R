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
parser$add_argument("--include", type = "character", nargs = "+", required = FALSE,
  help = "Include only specified feature directory.")
parser$add_argument("--exclude", type = "character", nargs = "+", required = FALSE,
  help = "Exclude specified feature directory.")
args <- parser$parse_args()
attach(args)

# original all features 
message("constructing loop and fet objects")
din <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/Spec_H3K27AC_ESC"
hichip <- file.path(dirname(din), "doc", "Spec_H3K27AC_ESC_LoopType.txt")
dout <- gsub("Spec_H3K27AC_ESC", "conFet", din)
exclude <- 'POLA2_TransTSS'
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = FALSE)
fet_fs <- dir(din, "*_heatmap.mat", recursive = TRUE, full.names = TRUE)
if (!is.null(include)) {
	includeStr <- paste0('/(', paste(include, collapse='|'), ')/')
	fet_fs <- fet_fs[grepl(includeStr, fet_fs)]
}
if (!is.null(exclude)) {
	excludeStr <- paste0('/(', paste(exclude, collapse='|'), ')/')
	fet_fs <- fet_fs[!grepl(excludeStr, fet_fs)]
}
fet.obj <- fetConst(fet_fs, small = 0.05) 
mViolinPlot(loop.obj, fet.obj, dout = dout)

# original all features 
message("constructing loop and fet objects")
din <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/Spec_H3K27AC_ESC"
hichip <- file.path(dirname(din), "doc", "Spec_H3K27AC_ESC_LoopType.txt")
dout <- gsub("Spec_H3K27AC_ESC", "conFetProm", din)
include <- 'POLA2_TransTSS'
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = FALSE)
fet_fs <- dir(din, "*_heatmap.mat", recursive = TRUE, full.names = TRUE)
if (!is.null(include)) {
	includeStr <- paste0('/(', paste(include, collapse='|'), ')/')
	fet_fs <- fet_fs[grepl(includeStr, fet_fs)]
}
if (!is.null(exclude)) {
	excludeStr <- paste0('/(', paste(exclude, collapse='|'), ')/')
	fet_fs <- fet_fs[!grepl(excludeStr, fet_fs)]
}
fet.obj <- fetConst(fet_fs, small = 0.05) 
mViolinPlotProm(loop.obj, fet.obj, dout = dout)