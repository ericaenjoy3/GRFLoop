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
  loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = ifelse(loopType|conHub, FALSE, TRUE), filterDist = FALSE)
}

if (length(hichip) >1) {
  loop.obj.list <- lapply(hichip, function(f)loopConst(f, score_col = NULL, filterUnknown = FALSE))
  cluster <- c()
  for (i in seq_along(loop.obj.list)) {
    string <- readline(prompt = paste0("please enter", i, "th cluster name: "))
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
  message("coreg")
  root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/coreg"
  info.obj <- infoConst()
  info.obj <- ProteinCodingInfo(info.obj)
  info.obj <- TPMInfo(info.obj)
  obj.list <- infoFilter(loop.obj, info.obj = info.obj)
  loop.obj <- obj.list[["loop.obj"]]
  message("linePlot")
  for (vtype in c("Prom", "Enh")) {
    linePlot(loop.obj, info.obj, 
      pdffout = file.path("/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/linePlot",
        gsub("PromEnh.txt", paste0(vtype, "_Line.pdf"), hichip)), 
      vtype = vtype, uniqueLoopGene = TRUE)
  }
  # info.obj <- geneCor(info.obj)
  # for (k in 2:4) {
  #   message("k ", k)
  #   shufPlot(loop.obj, info.obj, nmin = k , nmax =k, dout = root,
  #     tadStatpdf = file.path(root, paste0("Enh_nmin", k, "_nmax", k, "_tadStats.pdf")),
  #     coregBoxpdf = file.path(root, paste0("En_nmin", k, "_nmax", k, "_coregBox.pdf")),
  #     gcorBoxpdf = file.path(root, paste0("Enh_nmin", k, "_nmax", k, "_gcorBox.pdf")),
  #     glabBarpdf = file.path(root, paste0("Enh_nmin", k, "_nmax", k, "_glabBar.pdf")),
  #     uniqueLoopGene = TRUE)
  # }
  nmax <- maxCoreg(loop.obj, info.obj, coregfout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_coreg_gene.txt", basename(hichip))))
  message("maximum contacts in coreg: ", nmax)
  nmax <- as.numeric(readline("enter nmax: "))
  shufPlotMulti(loop.obj, info.obj, nmin = 2, nmax = nmax, dout = root, 
    pdffout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_Enh_glabBar.pdf", hichip)), 
    fout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_coreg_pval.txt", hichip)), 
    uniqueLoopGene = TRUE)
}


