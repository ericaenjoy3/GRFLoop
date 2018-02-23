#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export shufPlotMulti
setGeneric(name = "shufPlotMulti",
  def = function(loop.obj, info.obj, nmin, nmax, dout, glabBarpdf, uniqueLoopGene = FALSE){
    standardGeneric("shufPlotMulti")
  }
)

#' @rdname shufPlotMulti-methods
setMethod(f = "shufPlotMulti",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, nmin, nmax, dout, glabBarpdf, uniqueLoopGene) {

    dir.create(dout, showWarnings = FALSE, recursive = TRUE)

    dat_list <- list()

    for (j in seq(nmin, nmax)) {

      message("initalizing elements of dat_list")
      dat_list[[length(dat_list)+1]] <- list()

      message("retrieve genes for enhancer hubs")
      gene_list <- geneListProc(loop.obj, info.obj, nmin = j, nmax = j, type = "Enh", uniqueLoopGene = uniqueLoopGene)
    
      # coordinates for max interval
      message("coordinate for max interval")
      coord <- rbindlist(lapply(gene_list, function(gs, info.obj){
        idx <- chmatch(gs, info.obj@gene[["gene"]])
        stopifnot(all(!is.na(idx)))
        unique(info.obj@gene[idx][, list(chr = chr, start = min(start), end = max(end))])
      }, info.obj = info.obj))

      # global permutations of windows
      message("global permutations of windows")
      genep_list <- coordShulf(coord, info.obj, dout, nshulf = 20, nmin = j, nmax = j)

      # permutations within TADs
      message("permutations within TADs")      
      genet_list <- inTADShulf(gene_list, info.obj)

      message("assigning gene_list, genep_list and genet_list to dat_list elements")  
      dat_list[[length(dat_list)]][['gene_list']] <- gene_list
      dat_list[[length(dat_list)]][['genep_list']] <- genep_list
      dat_list[[length(dat_list)]][['genet_list']] <- genet_list

    }




    # pariwise lab plots
    if (!is.null(info.obj@gcor)) {
      gene_lab <- gene2pairwiseLab(gene_list, info.obj)
      genep_lab <- gene2pairwiseLab(genep_list, info.obj)
      genet_lab <- gene2pairwiseLab(genet_list, info.obj)
      dat <- rbindlist(list(data.table(type = "Genuine", gene_lab),
        data.table(type = "Global Random", genep_lab),
        data.table(type = "In-TAD Random", genet_lab)), use.names = FALSE)
      dat <- dat[, .N, by = .(type, gene_lab)]
      dat[, pct := round(100*N/sum(N), digits = 2), by = type][, gene_lab := factor(gene_lab)]
    }

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