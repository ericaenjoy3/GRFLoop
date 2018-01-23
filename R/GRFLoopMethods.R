#' @include GRFLoopClass.R GRFLoopGeneric.R
setMethod(f = "rmNonVarRNA",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj) {
    idx <- grep("RNA", fet.obj@hash[["sms"]], ignore.case = TRUE)[1]
    kpt.idx <- which(apply(fet.obj@dat_list[[idx]], 1, var) > 0)
    stopifnot(length(kpt.idx) < nrow(fet.obj@dat_list[[idx]]))
    # update dat_list slot of fet class
    # update loop slot of loop class
    # update g slot of loop class
    # update split slot of loop class
    if (length(kpt.idx) > 0) {
      fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[kpt.idx,])
      validObject(fet.obj)
      loop.obj@loop <- loop.obj@loop[kpt.idx,]
      loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
      loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% loop.obj@loop[["loop"]]))
      loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))
      if (!is.null(loop.obj@split)) {
        loop.obj@split <- factor(as.character(loop.obj@split)[kpt.idx], levels = unique(as.character(loop.obj@split)[kpt.idx]))
      }
      validObject(loop.obj)
    }
    # row scale RNA data
    invisible(sapply(grep("RNA", fet.obj@hash[["sms"]], ignore.case = TRUE), function(i){
      fet.obj@dat_list[[i]] <<- data.table(t(scale(t(fet.obj@dat_list[[i]]))))
    }))
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)

#' @rdname quantRm-methods
setMethod(f = "quantRm",
  signature = c("loop", "fet", "logical"),
  definition = function(loop.obj, fet.obj, dedup){
    score <- E(loop.obj@g)$score
    names(score) <- E(loop.obj@g)$loop
    if (!dedup) {
      score <- score[loop.obj@loop[["loop"]]]
    }
    thresh <- as.numeric(quantile(score, probs = c(0.25, 0.75)))
    bottom_loop <- E(loop.obj@g)$loop[E(loop.obj@g)$score <= thresh[1]]
    top_loop <- E(loop.obj@g)$loop[E(loop.obj@g)$score >= thresh[2]]
    # update g slot of loop object
    loop.obj@g <- delete.edges(loop.obj@g, which(!E(loop.obj@g)$loop %in% c(bottom_loop, top_loop)))
    loop.obj@g <- delete.vertices(loop.obj@g, which(igraph::degree(loop.obj@g)<1))
    # update loop slot of loop object
    kpt.idx <- which(loop.obj@loop[["loop"]] %in% c(bottom_loop, top_loop))
    loop.obj@loop <- loop.obj@loop[kpt.idx, ]
    loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
    split <- rep(0, length(loop.obj@loop[["loop"]]))
    split[loop.obj@loop[["loop"]] %in% top_loop] <- 1
    loop.obj@split <- factor(split, levels = unique(split))
    validObject(loop.obj)
    # update dat_list slot of fet object
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[kpt.idx,])
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)

#' @rdname orderLoop-methods
setMethod(f = "orderLoop",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, sm_nm, order_method){
    order_method <- match.arg(order_method)
    stopifnot(all(sm_nm %in% fet.obj@hash[["sms"]]))
    order2 <- function(..., decreasing = TRUE){ order(..., decreasing=decreasing) }
    idx <- which(grepl(sm_nm, fet.obj@hash[["sms"]], ignore.case = TRUE) & grepl("Enh", fet.obj@hash[["grps"]], ignore.case = TRUE))
    if (order_method == "last_column") {
      row_idx <- order(fet.obj@dat_list[[idx]][[ncol(fet.obj@dat_list[[idx]])]], decreasing = TRUE)
    }
    if (order_method == "diff") {
      row_idx <- order(fet.obj@dat_list[[idx]][[ncol(fet.obj@dat_list[[idx]])]] - fet.obj@dat_list[[idx]][[1]], decreasing = TRUE)
    }
    if (order_method == "last_to_first") {
      row_idx <- do.call(order2, fet.obj@dat_list[[idx]][,ncol(fet.obj@dat_list[[idx]]):1, with = FALSE])
    }
    if (order_method == "quant") {
      row_idx <- do.call(order2, data.table(quant = as.numeric(loop.obj@split), fet.obj@dat_list[[idx]][,ncol(fet.obj@dat_list[[idx]]):1, with = FALSE]))
    }
    # update split slot of loop object
    loop.obj@split <- factor(as.numeric(loop.obj@split)[row_idx],
      levels = unique(as.numeric(loop.obj@split)[row_idx]))
    # update loop slot of loop object
    loop.obj@loop <- loop.obj@loop[row_idx]
    loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
    validObject(loop.obj)
    # update dat_list slot of fet object
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)dat[row_idx])
    validObject(fet.obj)
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)

#' @rdname ceilFet-methods
setMethod(f = "ceilFet",
  signature = c("fet"),
  definition = function(fet.obj) {
    # update dat_list slot of fet object
    idx_list <- split(seq_len(nrow(fet.obj@hash)), fet.obj@hash[["sms"]])
    rng_list <- lapply(seq_along(idx_list), function(i, idx_list, fet.obj){
      idx_vec <- idx_list[[i]]
      if (!grepl("H3K27AC", names(idx_list)[i], ignore.case = TRUE)) {
        vec <- range(unlist(fet.obj@dat_list[idx_vec]))
      } else {
        vec <- quantile(unlist(fet.obj@dat_list[idx_vec]), probs = c(0, 0.95))
        names(vec) <- NULL
      }
      invisible(sapply(idx_vec, function(idx){
        col_nm <- colnames(fet.obj@dat_list[[idx]])
        fet.obj@dat_list[[idx]][, (col_nm) := lapply(.SD, function(x){x[x > vec[2] ] <- vec[2]; return(x)}), .SDcols = col_nm]}))
      return(rbind(vec))}, idx_list = idx_list, fet.obj = fet.obj)
    names(rng_list) <- names(idx_list)
    validObject(fet.obj)
    return(list(fet.obj = fet.obj, rng_list = rng_list))
  }
)

#' @rdname coordShulf-methods
setMethod(f = "coordShulf",
  signature = c("data.table", "info"),
  definition = function(coord, info.obj, dout, nshulf, ncon){
    stopifnot(file.exists(dout))
    genomef <- "/home/liuyiyua/athena/Gencode/mm10/sequence/chrNameLengthKnownChr.genome"
    exclf <- path.expand("~/athena/blacklist/mm10-blacklist.bed")
    nd = file.path(dout, "tmp")
    dir.create(nd, showWarnings = FALSE, recursive = TRUE)
    coord <- coord[mixedorder(coord[["chr"]], coord[["start"]], decreasing = FALSE)]
    f0 <- file.path(nd, paste0("b", 0, ".bed"))
    write.table(coord, file = f0, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    coordp <- rbindlist(lapply(1:nshulf, function(j, coord, f0, dout, genomef, exclf){
      fp <- file.path(nd, paste0("b", j, ".bed"))
      cmd <- paste("bedtools shuffle -seed", j, "-excl", exclf, "-i", f0,
        "-g", genomef, ">", fp)
      message(cmd)
      system(cmd)
      dat <- fread(fp, header = FALSE)
      stopifnot(nrow(dat) == nrow(coord))
      return(dat)
    }, coord = coord, f0 = f0, dout = dout, genomef = genomef, exclf = exclf))
    setnames(coordp, c("chr", "start", "end"))
    setkeyv(coordp, colnames(coordp))
    dic <- foverlaps(coordp, info.obj@gene, nomatch = 0, which = TRUE)
    genep_list <- split(info.obj@gene[dic[["yid"]]][["gene"]], dic$xid)
    for (j in seq_along(nrow(coordp))) {
      if (!j %in% names(gene_list)) {
        genep_list[[j]] <- NA
      }
    }
    unlink(nd, recursive=TRUE)
    idx <- sapply(genep_list, length) >= ncon
    stopifnot(sum(idx) > 0)
    genep_list <- genep_list[idx]
    return(genep_list)
  }
)

setMethod(f = "gene2direction",
  signature = c("list", "info"),
  definition = function(gene_list, info.obj){
    idx <- sapply(gene_list, function(x)is.character(x))
    message(sum(idx), " intervals overlap with genes")
    message(sum(!idx), " intervals do not overlap with genes")
    gene_list <- gene_list[idx]
    # deg labels for these genes
    col_nm <- colnames(info.obj@gene)[grep("^DEG", colnames(info.obj@gene))]
    deg_list <- lapply(gene_list, function(gs, info.obj, col_nm){
      idx <- chmatch(gs, info.obj@gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      deg_dat <- info.obj@gene[idx, col_nm, with = FALSE]
      stopifnot(nrow(deg_dat) == length(gs))
      return(deg_dat)
    }, info.obj = info.obj, col_nm = col_nm)
    up_list <- lapply(deg_list, function(dd){
      data.table(direction = "up", dd[, lapply(.SD, function(x){
        up <- 100 * sum(tolower(x) == "up", na.rm = TRUE)/length(x);
        dn <- 100 * sum(tolower(x) == "down", na.rm = TRUE)/length(x);
        return(up)
      }), .SDcols = colnames(dd)])
    })
    dn_list <- lapply(deg_list, function(dd){
      data.table(direction = "dn", dd[, lapply(.SD, function(x){
        dn <- 100 * sum(tolower(x) == "down", na.rm = TRUE)/length(x);
        return(dn)
      }), .SDcols = colnames(dd)])
    })
    stopifnot(length(up_list) == length(dn_list))
    deg_pct_list <- lapply(seq_along(up_list), function(j){
      rbindlist(list(up_list[[j]], dn_list[[j]]))
    })
    return(deg_pct_list)
  }
)
