#' @include GRFLoopClass.R GRFLoopGeneric.R
setMethod(f = "shufPlot",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, nmin, nmax, dout, tadStatpdf, coregBoxpdf) {
    # dedup loop in the loop slot of loop.obj
    kpt.idx <- !duplicated(loop.obj@loop[["loop"]])
    loop_hash <- loop.obj@loop[kpt.idx]
    setkeyv(loop_hash, "loop")
    # identify incident loop
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$type == "Enh"]
    ed <- incident_edges(loop.obj@g, ve)
    idx <- sapply(ed, length) >= nmin & sapply(ed, length) <= nmax
    # extract PromGene from loop slot of loop.obj
    gene_list <- lapply(ed[idx], function(es, loop_hash){
      lp <- gsub("|", "_", as_ids(es), fixed = TRUE)
      gs <- loop_hash[lp][["PromGene"]]
      stopifnot(all(!is.na(gs)))
      return(gs)
    }, loop_hash = loop_hash)
    # deg labels
    deg_pct_list <- gene2direction(gene_list, info.obj) # use gene2direction
    deg_pct <- melt(rbindlist(deg_pct_list), id.vars = "direction")
    # max interval for grouped gene
    span_vec <- sapply(gene_list, function(gs, info.obj){
      idx <- chmatch(gs, info.obj@gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      start <- min(info.obj@gene[idx][["start"]])
      end <- max(info.obj@gene[idx][["end"]])
      span <- end - start
      return(span)
    }, info.obj = info.obj)
    # coordinate for max interval
    coord <- rbindlist(lapply(gene_list, function(gs, info.obj){
      idx <- chmatch(gs, info.obj@gene[["gene"]])
      stopifnot(all(!is.na(idx)))
      unique(info.obj@gene[idx][, list(chr = chr, start = min(start), end = max(end))])
    }, info.obj = info.obj))
    # global permutations of windows
    genep_list <- coordShulf(coord, info.obj, dout, nshulf = 20, nmin = nmin, nmax = nmax)
    # deg labels for global permutations
    degp_pct_list <- gene2direction(genep_list, info.obj) # use gene2direction
    degp_pct <- melt(rbindlist(degp_pct_list), id.vars = "direction")
    # permutations within TADs
    genet_list <- inTADShulf(gene_list, info.obj)
    degt_pct_list <- gene2direction(genet_list, info.obj) # use gene2direction
    degt_pct <- melt(rbindlist(degt_pct_list), id.vars = "direction")
    # TADshipPlot
    TADshipPlot(list(gene_list, genep_list, genet_list), info.obj,
      nms = c("Genuine", "Global Random", "In-TAD Random"), pdffout = tadStatpdf)
    dat <- rbindlist(list(data.table(type = "Genuine", deg_pct),
      data.table(type = "Global Random", degp_pct),
      data.table(type = "In-TAD Random", degt_pct)), use.names = FALSE)
    # boxplot
    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(subset(dat, dat[["variable"]] == "DEG_MEF.ESC"), aes(x = direction, y = value, fill = type)) +
      geom_boxplot(position = position_dodge(width = 0.9)) +
      labs(x = "", y = "%") +
      theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"),
        legend.position = "top")
    ggsave(coregBoxpdf, p1)
  }
)
