libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}

root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP/coreg"
hichipfs <- c("ESC", "MEF", "CONSTANT")
hichipfs <- file.path("~/athena/Andreas_H3K27AC_HICHIP/doc", paste0("Spec_H3K27AC_", hichipfs, "_LoopType_PromEnh.txt"))

for (hichip in hichipfs) {
	cat("\n\n")
	message(hichip)
	loop.obj <- loopConst(hichip, score_col = NULL, filterUnknown = TRUE, filterDist = FALSE)
	info.obj <- infoConst()
	info.obj <- ProteinCodingInfo(info.obj)
	info.obj <- TPMInfo(info.obj)
	obj.list <- infoFilter(loop.obj, info.obj = info.obj)
	loop.obj <- obj.list[["loop.obj"]]
	nmax <- maxCoreg(loop.obj, info.obj, coregfout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_coreg_gene.txt", basename(hichip))))
	message("maximum contacts in coreg: ", nmax)
	nmax <- as.numeric(readline("enter nmax: "))
	shufPlotMulti(loop.obj, info.obj, nmin = 2, nmax = nmax, dout = root, 
		pdffout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_Enh_glabBar.pdf", basename(hichip))), 
		fout = file.path(root, gsub("(Spec_H3K27AC_[^_]+).+", "\\1_coreg_pval.txt", basename(hichip))), 
		uniqueLoopGene = TRUE)
}
