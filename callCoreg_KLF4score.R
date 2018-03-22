libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}

root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_KLF4_HICHIP"
dout <- file(root, "coreg_KLF4score")

hichip <- file.path("~/athena/Andreas_H3K27AC_HICHIP/doc", "Spec_H3K27AC_ESC_LoopType_PromEnh.txt")


cat("\n\n")
message(hichip)
loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = TRUE, filterDist = FALSE)
info.obj <- infoConst()
info.obj <- ProteinCodingInfo(info.obj)
info.obj <- TPMInfo(info.obj)
obj.list <- infoFilter(loop.obj, info.obj = info.obj)
loop.obj <- obj.list[["loop.obj"]]
# keep only co-regulation loop
obj.list <- selCoreg(loop.obj, info.obj = info.obj)
loop.obj <- obj.list[["loop.obj"]]

fs <- file.path(root, paste0(c("DAY3", "DAY6", "ESC"), "_KLF4_50K_H3K27AC_LoopType.txt"))

loop.obj.list <- list()
for (f in fs) {
	loop.obj.list[[length(loop.obj.list) + 1]] <- loopConst(f, score_col = 11, filterUnknown = TRUE, filterDist = FALSE)
}
