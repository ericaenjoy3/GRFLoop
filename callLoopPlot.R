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
parser$add_argument("--conHub", action="store_true",
  help = "make connectivity dot plots.")
args <- parser$parse_args()
attach(args)

message("constructing loop and fet objects")
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = ifelse(loopType, FALSE, TRUE))

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

if (conHub) {
  message("constructing info object")
  info.obj <- infoConst()
  message("filter info object by protein coding")  
  info.obj <- ProteinCodingInfo(info.obj)
  message("filter info object by TPM threshold")
  info.obj <- TPMInfo(info.obj)
  message("infoFilter")
  obj.list <- infoFilter(loop.obj, info.obj = info.obj)
  loop.obj <- obj.list[["loop.obj"]]
  hubPlot(loop.obj, pdffout, minSampling = FALSE, PromEnh = FALSE)
  hubPlot(loop.obj, gsub(".pdf", "_minSampled.pdf", pdffout), minSampling = TRUE, PromEnh = FALSE)  
}


