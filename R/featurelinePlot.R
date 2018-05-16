#' @include GRFLoopClass.R

#' @export featurelinePlot
setGeneric(name = "featurelinePlot",
  def = function(loop.obj, fet.obj, dout){
    standardGeneric("featurelinePlot")
  }
)

#' @rdname mViolinPlot-methods
setMethod(f = "featurelinePlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, dout){

    dir.create(dout, showWarnings = FALSE, recursive = TRUE)

    # dedup loop slot based on loop, loc1, loc2 and the first occurance of rowid
    dat <- copy(loop.obj@loop[, c("loop", "loc1", "loc2", "rowid"), with = FALSE])
    dat <- copy(dat[,.SD[1], by = c("loop", "loc1", "loc2")])

    # highly and lowly connected vertex
    g <- loop.obj@g
    vec <- igraph::degree(g, which(V(g)$vtype == "Enh"))

    thresh <- 1:3
    # match vertex to loc1 or loc2 and retrieve the first occurance of rowid
    # select relevant files and relevant rows
    dat_list <- list()
    for (j in seq_along(thresh)) {

      v <- names(vec)[vec == thresh[j]]

      loc1_idxv <- loc2_idxv <- c()

      for (i in seq_along(v)) {

        if (i %%100 ==0 ) {message(i)}

        loc1_idx <- copy(dat[loc1 %in% v[i], rowid])
        loc2_idx <- copy(dat[loc2 %in% v[i], rowid])

        if (length(loc1_idx) > 0) {

          loc1_idxv <- c(loc1_idxv, loc1_idx[1])

        } else if (length(loc2_idx) > 0) {

          loc2_idxv <- c(loc2_idxv, loc2_idx[1])

        } else {

          error("failed to find loc1 or loc2 match.")

        }       
      }

      chip_loc1_idx <- copy(fet.obj@hash[loc == "Loc1", which = TRUE])
      chip_loc2_idx <- copy(fet.obj@hash[loc == "Loc2", which = TRUE])

      dat_list[[length(dat_list) + 1]] <- lapply(seq_along(chip_loc1_idx), function(k){
          copy(rbind(fet.obj@dat_list[[chip_loc1_idx[k]]][loc1_idxv], fet.obj@dat_list[[chip_loc2_idx[k]]][loc2_idxv]))
        })

    }

    stopifnot(nrow(dat_list[[1]]) == sum(vec <= thresh[1] | vec >= thresh[2]))

    # merge low and high into one data.table
    ndat_list <- lapply(seq_along(dat_list[[1]]), function(j) {
      tmp_list <- lapply(1:3, function(k) {
        data.table(split = k, copy(dat_list[[k]][[j]][, ncol(dat_list[[k]][[j]]), with = FALSE]))
      })
      nsize <- min(sapply(tmp_list, nrow))
      set.seed(888)
      rbindlist(lapply(tmp_list, function(dd){
        dd[if(nrow(dd) > nsize){sample(1:nrow(dd), nsize)} else {1:nrow(dd)}, ]
      }), use.names = FALSE, fill = FALSE)
    })
    names(ndat_list) <- copy(fet.obj@hash[loc == "Loc1", chip])

    cmp <- data.table(combn(1:3, 2))

    plist <- lapply(seq_along(ndat_list), function(j, ndat_list, dout, cmp) {
      message("start plot ", j)
      pdffout <- file.path(dout, paste0(names(ndat_list)[j], "_medianViolin.pdf"))
      dd <- ndat_list[[j]]
      setnames(dd, gsub("\\-", "_",  colnames(dd)))
      p1 <- ggviolin(dd, x = "split", y = colnames(ndat_list[[j]])[-1], fill = "split", 
      add = "boxplot", add.params = list(fill = "white"), 
      xlab = "", ylab = "", 
      legend.title = "", title = names(ndat_list)[j]) +
      stat_compare_means(comparisons = cmp, method = "wilcox.test")
      ggsave(filename = pdffout, p1)
      message("complete plot ", j)
      return(p1)
    }, ndat_list = ndat_list, dout = dout, cmp = cmp)

  }
)