# hichip (*_LoopType_*.txt): 0-based coordinates

loopConst <- function(loop_f, score_col, filterUnknown = TRUE, filterDist = FALSE) {
  # g slot: edge: loop, etype, dist, score, cluster
  # g slot: vertex: name, vtype
  # loop slot: data.table of loop, gene1, gene2 and rowid
  if (is.character(loop_f)) { 
    dat <- fread(loop_f, header = TRUE, na.string = c("", "NA", "N/A"))
  } else if ("data.table" %in% class(loop_f)) {
    dat <- copy(loop_f)
  }
  if (!is.null(score_col)) {
    score_nm <- colnames(dat)[score_col]
    nscore_nm <- "score"
    setnames(dat, score_nm, nscore_nm)
  } else {
    nscore_nm <- NULL
  }
  cluster_nm <- if (any(grepl("cluster", colnames(dat)))) {
    "cluster"
  } else {
    NULL
  }
  # filter columns
  mCols <- c("loc1Chr", "loc1Start", "loc1End", "loc2Chr", "loc2Start", "loc2End", "gene1", "gene2", "loc1type", "loc2type")
  setnames(dat, colnames(dat)[1 : length(mCols)], mCols)

  dat <- dat[, c(mCols, nscore_nm, cluster_nm), with = FALSE]
  dat[, c("loc1") := paste0(loc1Chr, ":", loc1Start, "-", loc1End)]
  dat[, c("loc2") := paste0(loc2Chr, ":", loc2Start, "-", loc2End)]
  dat[, c("rowid") := seq_len(nrow(dat))]
  dat[, c("dist") := ifelse(loc1Chr == loc2Chr, floor(abs(loc2End-loc1End)), NA)]
  # fileter loops by distance of 1 MB
  if (filterDist) {
    dat <- dat[which(dist <= 1e6)]
  }  
  # filter columns
  dat[, c("loop", "etype") := list(paste0(loc1, "|", loc2), paste0(loc1type, "|", loc2type))]
  # filter 'Unknown' in loc1type or loc2type
  if (filterUnknown) {
    dat <- dat[loc1type != "Unknown" & loc2type != "Unknown"]
  }
  # filter columns
  e_dat <- unique(dat[, c("loc1", "loc2", "loop", "etype", "dist", nscore_nm, cluster_nm), with = FALSE])
  v_dat <- rbind(data.table(name = dat[["loc1"]], vtype = dat[["loc1type"]]),
    data.table(name = dat[["loc2"]], vtype = dat[["loc2type"]])) %>% unique()
  g <- graph_from_data_frame(e_dat, directed = FALSE, vertices = v_dat)
  loop_slot <- copy(dat[, c("loop", "loc1", "loc2", "gene1", "gene2", "rowid"), with = FALSE])
  loop.obj <- new("loop", g = g, loop = loop_slot)
  return(loop.obj)
}

# genef: 0-based coordinates
# tadf: 0-based coordinates

infoConst <- function(
  genef = path.expand("~/athena/Gencode/mm10/annotation/gencode.vM6.annotation.gene.bed"),
  fcf = path.expand("~/athena/RNA/RNA_seq/DF5154_2017_08_25/salmon/Daf_DiffAna_OrderFlip.xls"),
  tpmf = path.expand("~/athena/RNA/RNA_seq/DF5154_2017_08_25/salmon/TPM_rd_merge.txt"),
  tadf = path.expand("~/athena/HIC/HIC_seq/APP/TAD_mm10.bed"), 
  p_val = 0.01, fc_num = 1.5) {
    # gene slot
    stopifnot(all(file.exists(c(genef))))
    gene <- fread(genef, header = FALSE)
    setnames(gene, c("chr", "start", "end", "gene"))
    gene[, start := start + 1]
    setkeyv(gene, c("gene"))
    gene <- gene[which(chr %in% paste0("chr", c(1:19)))]
    gene[, c("type") := gsub("[^|]+\\|[^|]+\\|([^|]+)", "\\1", gene)]
    # add DEG columns to gene data.table
    fc <- fread(fcf, header = TRUE)
    fc[, c("gene") := paste(gid, gname, gtype, sep = "|")]
    col_nm <- copy(colnames(fc)[grep("^DEG", colnames(fc))])
    for (j in grep("^DEG", colnames(fc))) {
      set(fc, 1:nrow(fc), j, "NDiff")
      up_idx <- which(fc[[j-3]] > fc_num & fc[[j-2]] < p_val)
      dn_idx <- which(fc[[j-3]] < -fc_num & fc[[j-2]] < p_val)
      set(fc, up_idx, j, "Up")
      set(fc, dn_idx, j, "Down")    
    }
    ridx <- copy(fc[["gene"]] %in% gene[["gene"]])
    fc <- fc[ridx]
    for (col in col_nm) {
      nfc <- copy(fc[fc[[col]] != "NDiff"])
      gene[nfc[["gene"]], c(col) := nfc[[col]]]
      message(col)
      print(table(gene[[col]]))
      cat("\n")
    }
    setkeyv(gene, c("chr", "start", "end"))
    # add tadid column to gene data.table
    tad <- fread(tadf, header = FALSE)
    setnames(tad, c("chr", "start", "end"))
    tad[, start := start + 1L]
    setkeyv(tad, c("chr", "start", "end"))
    tad <- tad[chr %in% paste0("chr", c(1:19))]
    tad[, c("mstart", "mend") := list(as.integer(start - floor((start - shift(end, 1L, type="lag"))/2)),
      as.integer(end + ceiling((shift(start, 1L, type = "lead") - end)/2))), by = "chr"]
    tad[is.na(mstart), `:=`(mstart = start)]
    tad[is.na(mend), `:=`(mend = end)]
    tad[, `:=`(start = NULL, end = NULL)]
    setnames(tad, c("mstart", "mend"), c("start", "end"))
    ntad <- copy(tad)
    ntad[, start := start - 1]
    write.table(ntad, file = gsub(".bed", "_nogap.bed", tadf), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
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
    setkeyv(gene, "gene")
    # add TPM columns to gene data.table
    tpm <- SepTPMCnt(tpmf)$tpm.grp
    tpm <- data.table(gene = rownames(tpm), tpm, key = "gene")
    gid <- copy(gene[["gene"]])
    tpm <- copy(tpm[gid, nomatch = 0])
    setnames(tpm, colnames(tpm)[colnames(tpm) != "gene"], paste0("tpm_", colnames(tpm)[colnames(tpm) != "gene"]))
    stopifnot(identical(tpm[,gene], gene[,gene]))
    gene <- data.table(gene, tpm[, -grep("gene", colnames(tpm)), with = FALSE])
    setkeyv(gene, c("chr", "start", "end"))
    info.obj <- new("info", gene = gene, tad = tad)
    return(info.obj)
}

fetConst <- function(fet_fs, small = 0.05) {
  chip <- basename(dirname(fet_fs))
  loc <- gsub(".+(Loc\\d).+", "\\1", basename(fet_fs))
  hash <- data.table(chip = chip, loc = loc)
  hash[, (colnames(hash)) := lapply(.SD, function(vec){factor(vec, levels = unique(vec))}) , .SDcols = colnames(hash)]
  dat_list <- lapply(seq_along(fet_fs), function(i, fet_fs, chip, small){
    message("reading file ", fet_fs[i])
    dat <- fread(fet_fs[i], header = TRUE)
    if (grepl("RNA|PROMOTER_GENE", chip[i])) {
      dat[, c(colnames(dat)) := lapply(.SD, function(x){log2(x + small)}), .SDcols = colnames(dat)]
    }
    return(dat)
  }, fet_fs = fet_fs, chip = hash[["chip"]], small = small)
  names(dat_list) <- copy(hash[["chip"]])
  fet.obj <- new("fet", dat_list = dat_list, hash = hash)
  return(fet.obj)
}
