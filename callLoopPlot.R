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
parser$add_argument("--splithichip", action="store_true",
  help = "split hichip file based on cluster column.")
parser$add_argument("--pdffout", type = "character", required = FALSE,
  help = "output pdf")
parser$add_argument("--loopType", action="store_true",
  help = "make loop type bar plots.")
parser$add_argument("--loopDist", action="store_true",
  help = "make loop distance bar plots.")
parser$add_argument("--conHub", action="store_true",
  help = "make connectivity dot plots.")
parser$add_argument("--coreg", action="store_true",
  help = "make co-regulation plots.")
args <- parser$parse_args()
attach(args)

if (length(hichip) == 1 & !splithichip) {
  message("constructing loop and fet objects")
  loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = ifelse(loopType, FALSE, TRUE))
}

if (length(hichip) >1) {
  loop.obj.list <- lapply(hichip, function(f)loopConst(f, score_col = NULL, filterUnknown = FALSE))
  cluster <- c()
  for (i in seq_along(loop.obj.list)) {
    string <- readInput(prompt = paste0("please enter", i, "th cluster name: "))
    cluster <- c(cluster, string)
  }
  message("done with entering cluster names.")
  names(loop.obj.list) <- cluster
}

if (splithichip) {
  dat <- fread(hichip, header = TRUE)
  stopifnot("cluster" %in% colnames(dat))
  dat_list <- split(dat, dat[["cluster"]])
  loop.obj.list <- lapply(dat_list, function(dd)loopConst(dd, score_col = NULL, filterUnknown = FALSE))
}

if (loopType & length(hichip) == 1 & !splithichip) {
	message("loopTypePlot")
	loopTypePlot(loop.obj, pdffout)
	message("done with loopTypePlot")
}

if (loopDist & length(hichip) == 1 & !splithichip) {
  message("loopDistPlot")
  loopDistPlot(loop.obj, pdffout)
  message("done with loopDistPlot")	
}

if (conHub & length(hichip) > 1) {
  hubPlot(loop.obj.list, pdffout, minSampling = FALSE, subType = FALSE)
  hubPlot(loop.obj.list, gsub(".pdf", "_minSampled.pdf", pdffout), minSampling = TRUE, subType = FALSE)
  hubPlot(loop.obj.list, gsub(".pdf", "_subType.pdf", pdffout), minSampling = FALSE, subType = TRUE)
  hubPlot(loop.obj.list, gsub(".pdf", "_subType_minSampled.pdf", pdffout), minSampling = TRUE, subType = TRUE)  
}

if (coreg) {
  info.obj <- infoConst()
  info.obj <- ProteinCodingInfo(info.obj)
  info.obj <- TPMInfo(info.obj)
  obj.list <- infoFilter(loop.obj, fet.obj = NULL, info.obj = info.obj)
  loop.obj <- obj.list[["loop.obj"]]
  
}

