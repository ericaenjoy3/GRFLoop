#!/usr/bin/env Rscript

libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}

# OP either connect or heatmap
OP <- "connect"
type <- c("MEF", "ESC", "CONSTANT") #
root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP"
score_col <- NULL
info.obj <- infoConst()
info.obj <- ProteinCodingInfo(info.obj)
info.obj <- TPMInfo(info.obj)

for (i in seq_along(type)) {
  loop_f <- dir(file.path(root, "doc"), paste0("Spec_", type[i], ".txt"), recur = TRUE, full = TRUE)
  # constructing loop and fet objects
  message("constructing loop and fet objects")
  loop.obj <- loopConst(loop_f, score_col)
  message("infoFilter")
  obj.list <- infoFilter(loop.obj = loop.obj, info.obj = info.obj)
  loop.obj <- obj.list[["loop.obj"]]
  if (tolower(OP) == "connect") {
    conPlot(loop.obj = loop.obj, dout = file.path(root, paste0("Spec_", type[i])))
    for (k in 2) {
    	message("i ", i, " k ", k)
    	shufPlot(loop.obj, info.obj, nmin = k , nmax =k, dout = file.path(root, paste0("Spec_", type[i])),
    		tadStatpdf = file.path(file.path(root, paste0("Spec_", type[i])), paste0("Spec_", type[i], "_nmin", k, "_nmax", k, "_tadStats.pdf")),
    		coregBoxpdf = file.path(file.path(root, paste0("Spec_", type[i])), paste0("Spec_", type[i], "_nmin", k, "_nmax", k, "_coregBox.pdf")))
	  }
  }
}
