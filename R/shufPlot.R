#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export shufPlot
setGeneric(name = "shufPlot",
  def = function(loop.obj, info.obj, nmin, nmax, dout, tadStatpdf, coregBoxpdf, gcorBoxpdf, glabBarpdf, uniqueLoopGene = FALSE){
    standardGeneric("shufPlot")
  }
)

#' @rdname shufPlot-methods
setMethod(f = "shufPlot",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, nmin, nmax, dout, tadStatpdf, coregBoxpdf, gcorBoxpdf, glabBarpdf, uniqueLoopGene) {

    dir.create(dout, showWarnings = FALSE, recursive = TRUE)

    gene_list <- geneListProc(loop.obj, info.obj, nmin, nmax, type = "Enh", uniqueLoopGene = uniqueLoopGene)
    
    # deg labels
    deg_pct_list <- gene2direction(gene_list, info.obj) # use gene2direction
    deg_pct <- melt(rbindlist(deg_pct_list), id.vars = "direction")
    
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

    # boxplot of co-regulation labels
    theme_set(theme_grey(base_size = 15))
    if (nmin != nmax) {
      p1 <- ggplot(subset(dat, dat[["variable"]] == "DEG_MEF.ESC"), aes(x = direction, y = value, fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.9)) +
        labs(x = "", y = "%") +
        theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"),
          legend.position = "top")
    } else {
      dat[value <50, cate := "<50"]
      dat[value >= 50, cate := ">= 50"]
      dat[, value := NULL]
      dat <- dat[, N := .N, by = .(type, direction, variable, cate)] %>% unique()
      dat[, pct :=  100*N/sum(N), by = .(type, direction, variable)]
      p1 <- ggplot(subset(dat, dat[["variable"]] == "DEG_MEF.ESC"), aes(x = cate, y = pct, fill = type)) +
        geom_bar(stat="identity", position = position_dodge(width = 0.9)) +
        labs(x = "", y = "%") +
        theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"),
          legend.position = "top") + facet_grid(. ~ direction) + 
        geom_text(aes(label = N), hjust = 0.5, vjust = -0.5, size = 4, position = position_dodge(width = 0.9))
    }
    ggsave(coregBoxpdf, p1)

    # correlation coefficients
    if (!is.null(info.obj@gcor)) {
      gene_cor <- gene2pairwiseCor(gene_list, info.obj)
      genep_cor <- gene2pairwiseCor(genep_list, info.obj)
      genet_cor <- gene2pairwiseCor(genet_list, info.obj)
      dat <- rbindlist(list(data.table(type = "Genuine", gene_cor),
        data.table(type = "Global Random", genep_cor),
        data.table(type = "In-TAD Random", genet_cor)), use.names = FALSE)
      cmp <- data.table(combn(unique(dat[, type]), 2))
      p1 <- ggboxplot(dat, x = "type", y = "gene_cor", color = "type", palette = "jco", add = "jitter",
        xlab = "", ylab = "Spearman Correlation", legend.title = "") +
        stat_compare_means(comparison = cmp, method = "wilcox.test", label = "p.format")
      ggsave(gcorBoxpdf, p1)
    }

    # pariwise lab plots
    gene_lab <- gene2pairwiseLab(gene_list, info.obj)
    genep_lab <- gene2pairwiseLab(genep_list, info.obj)
    genet_lab <- gene2pairwiseLab(genet_list, info.obj)
    dat <- rbindlist(list(data.table(type = "Genuine", gene_lab),
      data.table(type = "Global Random", genep_lab),
      data.table(type = "In-TAD Random", genet_lab)), use.names = FALSE)
    dat <- dat[, .N, by = .(type, gene_lab)]
    dat[, pct := round(100*N/sum(N), digits = 2), by = type][, gene_lab := factor(gene_lab)]

    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(dat, aes(x = type, y = pct, fill = gene_lab))+
    geom_bar(stat = "identity") +
    labs(x = "",y = "%") +
    theme(legend.title = element_blank(),panel.spacing = unit(2, "lines"), legend.position = "top", 
      axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE))
    ggsave(glabBarpdf, p1)
  }
)

