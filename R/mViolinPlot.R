#' @include GRFLoopClass.R

#' @export mViolinViolin
setGeneric(name = "mViolinPlot",
  def = function(loop.obj, fet.obj, dout){
    standardGeneric("mViolinPlot")
  }
)

#' @rdname mViolinPlot-methods
setMethod(f = "mViolinPlot",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, dout){
    dir.create(dout, showWarnings = FALSE, recursive = TRUE)

    # dedup loop slot based on loop, loc1, loc2 and the first occurance of rowid
    dat <- loop.obj@loop[, c("loop", "loc1", "loc2", "rowid"), with = FALSE]
    dat <- dat[,.SD[1], by = c("loop", "loc1", "loc2")]

    # highly and lowly connected vertex
    g <- loop.obj@g
    vec <- igraph::degree(g)

    thresh <- c(1, 4)
    # match vertex to loc1 or loc2 and retrieve the first occurance of rowid
    # select relevant files and relevant rows
    dat_list <- list()
    for (j in seq_along(thresh)) {

      v <- if (j == 1) {
        names(vec)[vec <= thresh[j]]
      } else {
        names(vec)[vec >= thresh[j]]
      }

      loc1_idxv <- loc2_idxv <- c()
      for (i in seq_along(v)) {
        if (i %%100 ==0 ) {message(i)}
        loc1_idx <-dat[loc1 %in% v[i], rowid]
        loc2_idx <-dat[loc2 %in% v[i], rowid]
        if (length(loc1_idx) > 0) {
          loc1_idxv <- c(loc1_idxv, loc1_idx[1])
        } else if (length(loc2_idx) > 0) {
          loc2_idxv <- c(loc2_idxv, loc2_idx[1])
        } else {
          error("failed to find loc1 or loc2 match.")
        }       
      }

      chip_loc1_idx <- fet.obj@hash[loc == "Loc1", which = TRUE]
      chip_loc2_idx <- fet.obj@hash[loc == "Loc1", which = TRUE]

      dat_list[[length(dat_list) + 1]] <- lapply(seq_along(chip_loc1_idx), function(k){
          rbind(fet.obj@dat_list[[chip_loc1_idx[k]]][loc1_idxv], fet.obj@dat_list[[chip_loc2_idx[k]]][loc2_idxv])
        })

    }
    stopifnot(nrow(dat_list[[1]]) == sum(vec <= thresh[1] | vec >= thresh[2]))
    # merge low and high into one data.table
    ndat_list <- lapply(seq_along(dat_list[[1]]), function(j) {
      d1 <- data.table(split = "lowCon", dat_list[[1]][[j]][, ncol(dat_list[[1]][[j]]), with = FALSE])
      d2 <- data.table(split = "hiCon", dat_list[[2]][[j]][, ncol(dat_list[[2]][[j]]), with = FALSE])
      nsize <- min(nrow(d1), nrow(d2))
      set.seed(888)
      rbind(d1[if(nrow(d1) > nsize){sample(1:nrow(d1), nsize)} else {1:nrow(d1)}, ], 
        d2[if(nrow(d2) > nsize){sample(1:nrow(d2), nsize)} else {1:nrow(d2)}, ])
    })
    names(ndat_list) <- fet.obj@hash[loc == "Loc1", chip]

    cmp <- data.table(combn(c("lowCon", "hiCon"), 2))

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