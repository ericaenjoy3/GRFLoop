#!/usr/bin/env Rscript

libfs <- dir(path.expand("~/ComScripts/RPack/GRFLoop/R"), "*.R", recursive = TRUE, full.names = TRUE)
idx <- grepl("class", libfs, ignore.case = TRUE)
libfs <- libfs[c(which(idx), which(!idx))]
for (f in libfs) {
  message("source ", f)
  source(f)
}
# OP either connect or heatmap
OP <- "heatmap"
enhs <- c("SuperEnh", "NormalEnh") #
root <- "/athena/apostoloulab/scratch/liuyiyua/Andreas_H3K27AC_HICHIP"
ord <- c("RNA", "KLF4", "ATAC", "H3K27AC", "YY1")
order_method <- "quant"
score_col <- 10
info.obj <- infoConst()
info.obj <- ProteinCodingInfo(info.obj)
info.obj <- TPMInfo(info.obj)
for (i in seq_along(enhs)) {
  enhfs <- rev(dir(file.path(root, enhs[i]), "*Enh_heatmap.mat", recur = TRUE, full = TRUE))
  promfs <- rev(dir(file.path(root, enhs[i]), "*Prom_heatmap.mat", recur = TRUE, full = TRUE))
  if (OP == "heatmap") {
    enhfs <- enhfs[grepl("/(RNA|ATAC|KLF4|H3K27AC|YY1)/", enhfs, ignore.case = TRUE)]
    promfs <- promfs[grepl("/(RNA|ATAC|KLF4|H3K27AC|YY1)/", promfs, ignore.case = TRUE)]
    names(enhfs) <- basename(dirname(enhfs))
    names(promfs) <- basename(dirname(promfs))
    stopifnot(all(names(enhfs) %in% ord))
    stopifnot(all(names(promfs) %in% ord))
    enhfs <- enhfs[ord]
    promfs <- promfs[ord]
  }
  fet_fs <- c(enhfs, promfs)
  loop_f <- dir(file.path(root, "doc"), paste0("*", enhs[i], "_H3K27AC_ESC.txt"), recur = TRUE, full = TRUE)
  # constructing loop and fet objects
  message("constructing loop and fet objects")
  loop.obj <- loopConst(loop_f, score_col)
  fet.obj <- fetConst(fet_fs)
  if (tolower(OP) == "heatmap") {
    # quantile loops
    message("quantRm")
    obj.list <- quantRm(loop.obj, fet.obj, dedup = FALSE)
    loop.obj <- obj.list[["loop.obj"]]
    fet.obj <- obj.list[["fet.obj"]]
    # rm non variables looping
    message("rmNonVarRNA")
    obj.list <- rmNonVarRNA(loop.obj, fet.obj)
    loop.obj <- obj.list[["loop.obj"]]
    fet.obj <- obj.list[["fet.obj"]]
    # re-order loops
    message("orderLoop")
    obj.list <- orderLoop(loop.obj, fet.obj, sm_nm = "KLF4", order_method = order_method)
    loop.obj <- obj.list[["loop.obj"]]
    fet.obj <- obj.list[["fet.obj"]]
    # ceilling data values
    message("ceilling")
    dat_list <- ceilFet(fet.obj)
    # heatmap
    message("multiHeat")
    multiHeat(loop.obj, dat_list[["fet.obj"]], dat_list[["rng_list"]],
      pdffout = paste0(root, "/", enhs[i], "/", enhs[i], "_", order_method, "_heatmap.pdf"),
      scheme = ifelse(order_method == "last_column", "wr", "br"))
    # violin plot
    message("medianViolin")
    medianViolin(loop.obj, fet.obj,
      pdffout = paste0(root, "/", enhs[i], "/", enhs[i], "_", order_method, "_medianViolin.pdf"))
  }
  if (tolower(OP) == "connect") {
    message("infoFilter")
    obj.list <- infoFilter(loop.obj, fet.obj, info.obj)
    loop.obj <- obj.list[["loop.obj"]]
    fet.obj <- obj.list[["fet.obj"]]
    # conPlot(loop.obj, fet.obj, dout = root)
    shufPlot(loop.obj, info.obj, ncon = 4 , dout = root,
      pdffout = paste0(root, "/", enhs[i], "/", enhs[i], "_nmin", "_nmax", "_coregBox.pdf"))
  }
}
