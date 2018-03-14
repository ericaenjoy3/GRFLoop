#' @include GRFLoopClass.R

#' @export coordShulf
setGeneric(name = "coordShulf",
  def = function(coord, info.obj, dout, nshulf, nmin, nmax){
    standardGeneric("coordShulf")
  }
)

#' @rdname coordShulf-methods
setMethod(f = "coordShulf",
  signature = c("data.table", "info"),
  definition = function(coord, info.obj, dout, nshulf, nmin, nmax){
    stopifnot(file.exists(dout))
    genomef <- "/home/liuyiyua/athena/Gencode/mm10/sequence/autosome.genome"
    exclf <- path.expand("~/athena/blacklist/mm10-blacklist.bed")
    nd <- file.path(dout, "tmp")
    dir.create(nd, showWarnings = FALSE, recursive = TRUE)
    coord <- coord[mixedorder(coord[["chr"]], coord[["start"]], decreasing = FALSE)]
    f0 <- file.path(nd, paste0("b", 0, ".bed"))
    ncoord <- copy(coord)
    set(ncoord, i = 1:nrow(ncoord), j = 2, value = as.numeric(ncoord[[2]] - 1))
    write.table(coord, file = f0, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    coordp <- rbindlist(lapply(1:nshulf, function(j, coord, f0, dout, genomef, exclf){
      fp <- file.path(nd, paste0("b", j, ".bed"))
      cmd <- paste("bedtools shuffle -seed", j, "-excl", exclf, "-i", f0,
        "-g", genomef, ">", fp)
      message(cmd)
      system(cmd)
      dat <- fread(fp, header = FALSE)
      set(dat, i = 1:nrow(dat), j = 2, value = as.numeric(dat[[2]] + 1))
      stopifnot(nrow(dat) == nrow(coord))
      return(dat)
    }, coord = coord, f0 = f0, dout = dout, genomef = genomef, exclf = exclf))
    setnames(coordp, c("chr", "start", "end"))
    setkeyv(coordp, colnames(coordp))
    dic <- foverlaps(coordp, info.obj@gene, nomatch = 0, which = TRUE)
    genep_list <- split(info.obj@gene[dic[["yid"]]][["gene"]], dic$xid)
    for (j in seq_along(nrow(coordp))) {
      if (!j %in% names(genep_list)) {
        genep_list[[j]] <- NA
      }
    }
    unlink(nd, recursive = TRUE)
    idx <- sapply(genep_list, length) >= nmin & sapply(genep_list, length) <= nmax
    stopifnot(sum(idx) > 0)
    genep_list <- copy(genep_list[idx])
    return(genep_list)
  }
)
