#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export maxCoreg
setGeneric(name = "maxCoreg",
  def = function(loop.obj, info.obj, coregfout){
    standardGeneric("maxCoreg")
  }
)

#' @rdname geneList-methods
setMethod(f = "maxCoreg",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, coregfout) {

    type <- "Enh"

    # dedup loop in the loop slot of loop.obj
    kpt.idx <- copy(!duplicated(loop.obj@loop[["loop"]]))
    loop_hash <- copy(loop.obj@loop[kpt.idx])
    setkeyv(loop_hash, "loop")

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    
    # output gene coregulated by the same enhancer hub
    # columns:
    # 1, enh hub coordinates
    # 2, loop number
    # 3, promoter fragements,
    # 3, gene numbers
    # 4, unique gene numbers
    # 5, up-regulated unique gene numbers
    # 6, up-regulated unique gene pair numbers
    # 7, up-regulated unique genes
    # filter by unique gene numbers > 1
    dd <- rbindlist(lapply(seq_along(ed), function(i) {
      gene <- unlist(lapply(as_ids(ed[[i]]), function(string)loop.obj@loop[loop %in% string, gene1]))
      uniq_gene <- unique(gene)
      uniq_gene_stats <- data.table(gene)
      uniq_gene_stats <- paste(uniq_gene_stats[,.N, by = gene][, paste0(gene, "(", N, ")")], collapse = ",")
      deg_lab <- sapply(uniq_gene, function(string)info.obj@gene[gene %in% string, DEG_ESC.MEF])
      idx <- grepl("Up", deg_lab, ignore.case = TRUE)
      upreg_uniq_gene <- ifelse(sum(idx) > 0, paste(uniq_gene[idx], collapse = ","), "NA")
      upreg_uniq_gpair_num <- ifelse(sum(idx) >= 2, choose(sum(idx), 2), 0)
      tmp_dat <- data.table(enh = names(ed)[i],
        loop_num = length(ed[[i]]),
        prom = paste(gsub(paste0("\\|{0,1}", names(ed)[i], "\\|{0,1}"), "", as_ids(ed[[i]])), collapse = ","),
        gene_num = length(gene),
        uniq_gene_num = length(unique(gene)),
        upreg_uniq_gene_num = sum(idx),
        upreg_uniq_gpair_num = upreg_uniq_gpair_num,
        uniq_gene = uniq_gene_stats,
        upreg_uniq_gene = upreg_uniq_gene)
      return(tmp_dat)
    }))
    write.table(dd, file = coregfout, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

    # max connection for which the number of hub > 1
    maxco <- data.table(connum = sapply(ed, length))
    maxco <- maxco[, .N, by = connum]

    return(maxco[N>1, max(connum)])
  }
)