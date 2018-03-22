libfs <- dir(path.expand("~/athena/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (j in seq_along(libfs)) {
  message("source file ", j, " ",libfs[j])
  source(libfs[j])
}

root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_KLF4_HICHIP"
dout <- file.path(root, "coreg_KLF4score")
dir.create(dout, showWarnings = FALSE, recursive = TRUE)

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

sortlabs <- function(string) {
	vec <- c(gsub("(.+)\\|(.+)", "\\1", string), gsub("(.+)\\|(.+)", "\\2", string))
	paste(mixedsort(vec), collapse = "|")
}

dat_list[[1]] <- data.table(loop = sapply(as_ids(E(loop.obj@g)), sortlabs), key = "loop")

fs <- file.path(root, "doc", paste0(c("DAY3", "DAY6", "ESC"), "_KLF4_50K_H3K27AC_LoopType.txt"))

counter <- 1
for (f in fs) {
	counter <- counter + 1
	tmp.loop.obj <- loopConst(f, score_col = 11, filterUnknown = FALSE, filterDist = FALSE)
	ng <- intersection(tmp.loop.obj@g, loop.obj@g, keep.all.vertices = FALSE)
	dat_list[[counter]] <- data.table(loop = sapply(as_ids(E(ng)), sortlabs), score = E(ng)$score, key = "loop")
}
dat <- Reduce(function(...) merge(..., all = TRUE), dat_list)
dat[is.na(dat)] <- 0
setnames(dat, colnames(dat)[-1], c("DAY3", "DAY6", "ESC"))
ndat <- copy(dat)
ndat[ndat>0] <- 1

order2 <- function(..., decreasing = TRUE){ order(..., decreasing = decreasing) }
row_idx <- do.call(order2, ndat[,ncol(ndat):2, with = FALSE])

# bks <- c(0, seq(0, 10, 1), seq(20, dat[, max(.SD), .SDcols = colnames(dat)[-1]], 10))
bks <- c(0, 1)
colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Reds"))(length(bks))

ht_list <- ComplexHeatmap::Heatmap(dat[row_idx, !"loop", with = FALSE],
        col = circlize::colorRamp2(bks, colors),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        use_raster = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        combined_name_fun = NULL,
        name = "score",
        show_heatmap_legend = FALSE)

pdffout <- file.path(dout, "H3K27AC_Coreg_Loop_KLF4_Loop_Score.pdf")
pdf(pdffout)
draw(ht_list)
dev.off()
