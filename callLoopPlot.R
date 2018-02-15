#!/usr/bin/env Rscript

libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--hichip", type = "character", required = FALSE,
  help = "A query loop file with columns being locus 1 and 2, genes 1 and 2 and optional columns.")
parser$add_argument("--pdffout", type = "character", required = FALSE,
  help = "output pdf")
parser$add_argument("--loopType", action="store_true",
  help = "make loop type bar plots.")
parser$add_argument("--loopDist", action="store_true",
  help = "make loop distance bar plots.")
args <- parser$parse_args()
attach(args)

message("constructing loop and fet objects")
loop.obj <- loopConst(hichip, score_col = NULL)

if (loopType) {
	message("loopTypePlot")
	loopTypePlot(loop.obj, pdffout)
	message("done with loopTypePlot")
}

if (loopDist) {
  message("loopDistPlot")
  loopDistPlot(loop.obj, pdffout)
  message("done with loopDistPlot")	
}

