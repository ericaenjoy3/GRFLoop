#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export shufPlotMulti
setGeneric(name = "shufPlotMulti",
  def = function(loop.obj, info.obj, nmin, nmax, dout, pdffout, fout, uniqueLoopGene = FALSE){
    standardGeneric("shufPlotMulti")
  }
)

#' @rdname shufPlotMulti-methods
setMethod(f = "shufPlotMulti",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, nmin, nmax, dout, pdffout, fout, uniqueLoopGene) {

    dir.create(dout, showWarnings = FALSE, recursive = TRUE)

    dat_list <- list()
    dat_list[['gene_list']] <- list()
    dat_list[['genep_list']] <- list()
    dat_list[['genet_list']] <- list()

    # retrieve genes contacted by the same enhancer hub
    for (j in seq(nmin, nmax)) {

      message("initalizing elements of dat_list ", j, " th iteration")

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

      message("assigning gene_list, genep_list and genet_list to dat_list elements\n")  
      dat_list[['gene_list']][[length(dat_list[['gene_list']]) + 1]] <- gene_list
      dat_list[['genep_list']][[length(dat_list[['genep_list']]) + 1]] <- genep_list
      dat_list[['genet_list']][[length(dat_list[['genet_list']]) + 1 ]] <- genet_list
      for (s in c('gene_list', 'genep_list', 'genet_list')) {
        names(dat_list[[s]])[length(dat_list[[s]])] <- j
      }
    }

    # gene list to gene pairs
    # unpack gene_list, genep_list and genet_list
    # unpack enhancer hub classes (contacts 2, 3, ...)
    # unpack genes regulated by the same enhancer.
    gpair_list <- lapply(dat_list, function(num_list){
      dd_list <- lapply(num_list, function(g_list) {
        dd <- rbindlist(lapply(g_list, function(gs) {
          data.table(t(combn(gs, 2)))
        }))
        return(dd)
      })
      nd <-rbindlist(dd_list) %>% unique()
      return(nd)
    })

    glab_list <- gene2pairwiseLab(gpair_list, info.obj)
    glab <- rbindlist(lapply(seq_along(glab_list), function(i){
      type <- switch(names(glab_list)[i], gene_list = "Genuine", genep_list = "Global random", genet_list = "TAD-matched random", "NA")
      dd <- data.table(type = type, gene_lab = glab_list[[i]])
      dd <- dd[, .N, by = .(gene_lab, type)]
      dd[, pct := round(100*N/sum(N), digits = 2), by = .(type)]
      return(dd)
    }))
    glab <- glab[gene_lab !=0]
    glab[, pct := round(100*N/sum(N), digits = 2), by = .(type)]

    con <- file(fout, "w")
    idx_mat <- combn(glab[, unique(type)], 2)
    dd <- dcast(glab[, c("type", "gene_lab", "N")], type ~ gene_lab, value.var = "N")
    dd[is.na(dd)] <- 0
    write.table(dd, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", file = con, append = TRUE)
    cat("\n\n", sep = "", file = con)
    for (j in 1:ncol(idx_mat)) {
      p_val <- fisher.test(as.matrix(dd[type %in% idx_mat[, j], !"type"]))$p.value
      cat("comparing ", dd[type %in% idx_mat[,j], paste(type, collapse = " vs ")], " p value: ", p_val, "\n\n", sep = "", file = con)
    }
    close(con)

    ndat<- glab
    ndat[, gene_lab := factor(gene_lab, levels = sort(unique(gene_lab)))]

    theme_set(theme_grey(base_size=15))
    p1 <- ggplot(ndat, aes(type, pct, fill = gene_lab)) + 
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    labs(x = "",y = "%") +
    theme(legend.title = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "top", 
      axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE))
    
    ggsave(pdffout, p1)

  }
)