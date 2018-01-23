loopConst <- function(loop_f, score_col) {
  dat <- fread(loop_f, header = TRUE)
  score_nm <- colnames(dat)[score_col]
  dat <- dat[, c("PromChr", "PromStart", "PromEnd", "EnhChr", "EnhStart", "EnhEnd", "PromGene", score_nm), with = FALSE]
  setnames(dat, score_nm, "score")
  dat[, c("Prom") := paste0(PromChr, ":", PromStart, "-", PromEnd)]
  dat[, c("Enh") := paste0(EnhChr, ":", EnhStart, "-", EnhEnd)]
  dat[, c("rowid") := seq_len(nrow(dat))]
  dat[, c("PromChr", "PromStart", "PromEnd", "EnhChr", "EnhStart", "EnhEnd") := NULL]
  dat <- dat[, c("Prom", "Enh", "PromGene", "score", "rowid"), with = FALSE]
  dat[, c("loop") := paste0(Prom, "_", Enh)]
  e_dat <- unique(dat[, c("Prom", "Enh", "loop", "score"), with = FALSE])
  v_dat <- unique(rbind(data.frame(V = dat$Prom, type = "Prom"),
    data.frame(V = dat$Enh, type = "Enh")))
  g <- graph_from_data_frame(e_dat, directed = FALSE, vertices = v_dat)
  loop_slot <- dat[, c("loop", "PromGene", "rowid"), with = FALSE]
  loop.obj <- new("loop", g = g, loop = loop_slot)
  return(loop.obj)
}

infoConst <- function(
  genef = path.expand("~/athena/Gencode/mm10/annotation/gencode.vM6.annotation.gene.bed"),
  fcf = path.expand("~/athena/RNA/RNA_seq/DF5154_2017_08_25/hera/Daf_DiffAna_OrderFlip.xls"),
  tadf = path.expand("~/athena/HIC/HIC_seq/APP/TAD_mm10.bed")) {
    # gene slot
    stopifnot(all(file.exists(c(genef))))
    gene <- fread(genef, header = FALSE)
    setnames(gene, c("chr", "start", "end", "gene"))
    setkeyv(gene, c("chr", "start", "end"))
    gene[, c("type") := gsub("[^|]+\\|[^|]+\\|([^|]+)", "\\1", gene)]
    fc <- fread(fcf, header = TRUE)
    fc[, c("gene") := paste(gid, gname, gtype, sep = "|")]
    col_nm <- colnames(fc)[grep("^DEG", colnames(fc))]
    for (col in col_nm) {
      nfc <- fc[fc[[col]] != "NDiff"]
      idx <- chmatch(nfc[["gene"]], gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      gene[idx, c(col) := nfc[[col]]]
      message(col)
      print(table(gene[[col]]))
      cat("\n")
    }
    tad <- fread(tadf, header = FALSE)
    setnames(tad, c("chr", "start", "end"))
    setkeyv(tad, c("chr", "start", "end"))
    tad[, c("mstart", "mend") := list(as.integer(start - floor((start - shift(end, 1L, type="lag"))/2)),
      as.integer(end + ceiling((shift(start, 1L, type = "lead") - end)/2))), by = "chr"]
    tad[is.na(mstart), `:=`(mstart = start)]
    tad[is.na(mend), `:=`(mend = end)]
    tad[, `:=`(start = NULL, end = NULL)]
    setnames(tad, c("mstart", "mend"), c("start", "end"))
    fwrite(tad, file = gsub(".bed", "_nogap.bed", tadf), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    setkeyv(tad, c("chr", "start", "end"))
    # filter gene and tad based on common chr
    tad <- tad[chr %in% gene[["chr"]]]
    gene <- gene[chr %in% tad[["chr"]]]
    tad[, c("tadid") := seq_len(nrow(tad))]
    gene[, c("tadid") := as.integer(NA)]
    dic <- foverlaps(gene, tad, nomatch = 0, which = TRUE)
    dic <- unique(dic[, list(nxid = unique(xid), nyid = yid[1]),  by = "xid"][, c("nxid", "nyid"), with = FALSE])
    setnames(dic, c("nxid", "nyid"), c("xid", "yid"))
    gene[dic$xid, c("tadid") := tad[["tadid"]][dic$yid]]
    info.obj <- new("info", gene = gene)
}

fetConst <- function(fet_fs, small = 0.05) {
  hash <- data.table(
    sms = gsub("^[^_]+_([^_]+)_.*", "\\1", basename(fet_fs)),
    grps = gsub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", basename(fet_fs)),
    enhs = basename(dirname(dirname(fet_fs)))
  )
  hash[, c("grps_enhs") := paste(grps, enhs, sep = "_")]
  lapply(seq_len(ncol(hash)), function(i){
    hash[[i]] <<- factor(hash[[i]], levels = unique(hash[[i]]));
    invisible(NULL)
  })
  dat_list <- lapply(seq_along(fet_fs), function(i, fet_fs, sms, small){
    message("reading file ", fet_fs[i])
    dat <- fread(fet_fs[i], header = TRUE)
    if (grepl("RNA", sms[i])) {
      dat[, c(colnames(dat)) := lapply(.SD, function(x){log2(x + small)}), .SDcols = colnames(dat)]
    }
    return(dat)
  }, fet_fs = fet_fs, sms = hash[["sms"]], small = small)
  names(dat_list) <- hash[["sms"]]
  fet.obj <- new("fet", dat_list = dat_list, hash = hash)
  return(fet.obj)
}
