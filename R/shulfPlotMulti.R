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

    Tpair_list <- list()
    invisible(sapply(1:3, function(j)Tpair_list[[j]] <<- list()))
    # retrieve genes contacted by the same enhancer hub
    for (j in seq(nmin, nmax)) {

      message("initalizing elements of dat_list ", j, " th iteration")

      message("genePair")
      gpair_list <- genePair(loop.obj, info.obj, type = "Enh", nmin = j, nmax = j)
      for (i in seq(gpair_list)) {
      	if (length(Tpair_list[[i]]) == 0) {
      		Tpair_list[[i]] <- gpair_list[[i]]
      	} else {
        	Tpair_list[[i]] <- rbind(Tpair_list[[i]], gpair_list[[i]])
        }
      }
    }
    sapply(seq(Tpair_list), function(j)Tpair_list[[j]] <<- unique(Tpair_list[[j]]))

    # gene pair to DEG labels
    Tlab_list <- list()
    for (j in seq(Tpair_list)) {
    	Tlab_list[[j]] <- gene2pairwiseLab(Tpair_list[[j]], info.obj)
    }

    glab <- rbindlist(lapply(seq_along(Tlab_list), function(i){
      type <- switch(i, '1' = "Genuine", '2' = "Global random", '3' = "TAD-matched random", "NA")
      dd <- data.table(type = type, gene_lab = Tlab_list[[i]])
      dd <- dd[, .N, by = .(gene_lab, type)]
      dd[, pct := round(100*N/sum(N), digits = 2), by = .(type)]
      return(copy(dd))
    }))

    con <- file(fout, "w")

    dd <- copy(dcast(glab[, c("type", "gene_lab", "N")], type ~ gene_lab, value.var = "N"))
    dd[is.na(dd)] <- 0

    write.table(dd, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", file = con, append = TRUE)
    cat("\n\n", sep = "", file = con)

    ndat <- copy(glab[gene_lab !=0])
    idx_mat <- combn(ndat[, unique(type)], 2)	
    ndat[, pct := round(100*N/sum(N), digits = 2), by = .(type)]
    dd <- copy(dd[, !'0', with = FALSE])

    for (j in 1:ncol(idx_mat)) {

      p_val <- fisher.test(as.matrix(dd[type %in% idx_mat[, j], !"type"]))$p.value
      cat("comparing ", dd[type %in% idx_mat[,j], paste(type, collapse = " vs ")], " p value: ", p_val, "\n\n", sep = "", file = con)

    }
    close(con)

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