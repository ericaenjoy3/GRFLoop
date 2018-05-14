#!/usr/bin/env Rscript

libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}
library(argparse)

reportFunc <- function(loop.obj) {
	type = "Enh"
	stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
	kpt_num <- sapply(ed, length) >= 2
	message(sum(kpt_num), " enhancer fragments contacting two or more promoter fragments.")
}

f <- "Spec_H3K27AC_ESC_LoopType_PromEnh.txt"
hichip <- file.path("/home/liuyiyua/athena/Andreas_H3K27AC_HICHIP/doc", f)
message("constructing loop and fet objects")
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = TRUE, filterDist = FALSE)
reportFunc(loop.obj)
