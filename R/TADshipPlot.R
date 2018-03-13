#' @include GRFLoopClass.R

#' @export TADshipPlot
setGeneric(name = "TADshipPlot",
  def = function(gene_list_list, info.obj, nms, pdffout){
    standardGeneric("TADshipPlot")
  }
)

#' @rdname TADshipPlot-methods
setMethod(f = "TADshipPlot",
  signature = c("list", "info"),
  definition = function(gene_list_list, info.obj, nms, pdffout){
    stopifnot(length(nms) == length(gene_list_list))
    # barplot
    tadnum_list <- lapply(gene_list_list, function(gene_list){sapply(gene_list, function(gid){
        length(unique(info.obj@gene[gene %in% gid, tadid]))
      })
    })
    tadnum <- rbindlist(lapply(seq_along(tadnum_list), function(j)data.table(nms = nms[j], tadnum = tadnum_list[[j]])))
    tadnum <- tadnum[, .N, by = .(nms, tadnum)]
    tadnum <- tadnum[order(nms, -tadnum)]
    #
    tadhom_list <- lapply(gene_list_list, function(gene_list){sapply(
      gene_list, function(gid) {
        tadid <- copy(info.obj@gene[gene %in% gid, tadid])
        freq_tadid <- as.numeric(tail(names(sort(table(tadid))), 1))
        tadhom <- 100 * sum(tadid == freq_tadid)/length(tadid)
        return(tadhom)
      })
    })
    tadhom <- rbindlist(lapply(seq_along(tadhom_list), function(j)data.table(nms = nms[j], tadhom = tadhom_list[[j]])))
    p1 <- ggbarplot(tadnum, x = "tadnum", y = "N", fill = "nms", color = "white", palette = "jco", rotate = TRUE,
      position = position_dodge(0.9), xlab = "TADs per Hub", ylab = "Counts", legend.title = "", ggtheme = theme_minimal())
    cmp <- data.table(combn(nms, 2))
    p2 <- ggboxplot(tadhom, x = "nms", y = "tadhom", color = "nms", palette = "jco", add = "jitter",
      xlab = "", ylab = "TAD Homogeneity per Hub (%)", legend.title = "") +
      stat_compare_means(comparison = cmp, method = "wilcox.test", label = "p.format")
    pdf(pdffout, width = 7*2)
    multiplot(p1, p2, cols=2)
    dev.off()
  }
)
