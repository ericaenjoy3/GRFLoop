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

      message("assigning gene_list, genep_list and genet_list to dat_list elements\n")  
      dat_list[[length(dat_list)]][['gene_list']] <- gene_list
      dat_list[[length(dat_list)]][['genep_list']] <- genep_list
      dat_list[[length(dat_list)]][['genet_list']] <- genet_list
      names(dat_list)[length(dat_list)] <- j

    }

    glab_list <- lapply(nmin:nmax, function(j){
      message("processing ", j)
      message("retrieve gene list")
      idx <- which(names(dat_list) == j)
      stopifnot(length(idx) == 1)
      gene_list <- dat_list[[idx]][['gene_list']]
      genep_list <- dat_list[[idx]][['genep_list']]
      genet_list <- dat_list[[idx]][['genet_list']]
      message("starting gene2pairwiseLab")
      gene_lab <- gene2pairwiseLab(gene_list, info.obj)
      genep_lab <- gene2pairwiseLab(genep_list, info.obj)
      genet_lab <- gene2pairwiseLab(genet_list, info.obj)
      message("done with gene2pairwiseLab")
      message("making data.table")
      dat <- rbindlist(list(data.table(type = "Genuine", gene_lab),
        data.table(type = "Global Random", genep_lab),
        data.table(type = "In-TAD Random", genet_lab)), use.names = FALSE)
      dat <- dat[, .N, by = .(type, gene_lab)]
      dat[, pct := round(100*N/sum(N), digits = 2), by = type][, gene_lab := factor(gene_lab)]
      dat[, connum := j]
      message("done with data.table\n\n")
      return(dat)
    })

    glab <- rbindlist(glab_list, use.names = FALSE)
    glab <- glab[gene_lab != '0']
    con <- file(fout, "w")
    idx_mat <- combn(glab[, unique(type)], 2)
    for (num in glab[, unique(connum)]) {
      dd <- dcast(glab[connum==num, c("type", "gene_lab", "N")], type ~ gene_lab, value.var = "N")
      write.table(dd, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", file = con, append = TRUE)
      for (j in 1:ncol(idx_mat)) {
        p_val <- fisher.test(as.matrix(dd[type %in% idx_mat[, j], !"type"]))$p.value
      }
      cat("p value: ", p_val, "\n\n", sep = "", file = con)
    }
    close(con)

    p1 <- ggerrorplot(glab, x = "type", y = "pct", color = "gene_lab", palette = "jco", xlab = "", 
      ylab = "%", legend.title = "", add = "mean_se") +
      rotate_x_text(45)
    ggsave(pdffout, p1)

  }
)