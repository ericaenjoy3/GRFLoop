#' @include GRFLoopClass.R GRFLoopGeneric.R

#' @export genePair
setGeneric(name = "genePair",
  def = function(loop.obj, info.obj, dout, type, nmin, nmax, random.it = TRUE){
    standardGeneric("genePair")
  }
)

#' @rdname genePair-methods
setMethod(f = "genePair",
  signature = c("loop", "info"),
  definition = function(loop.obj, info.obj, dout, type, nmin, nmax, random.it) {

    # identify incident loop
    stopifnot(sum(V(loop.obj@g)$vtype == type) > 0)
    ve <- V(loop.obj@g)$name[V(loop.obj@g)$vtype == type]
    ed <- incident_edges(loop.obj@g, ve)
    kpt_num <- sapply(ed, length) >= nmin & sapply(ed, length) <= nmax

    # (genunine) extract PromGene from loop slot of loop.obj
    gene_list <- list()
    gene_list[[1]] <- lapply(ed[which(kpt_num)], function(es, loop_hash){
      # potentially multiple loops => mutiple fragments
      lp <- as_ids(es)
      dd <- loop_hash[loop %in% lp]
      stopifnot(all(dd[, all(!is.na(gene1)), by = loop][, V1]))
      dd <- dd[, unique(gene1), by = loop]
      gs <- split(dd[, V1], dd[, loop])
      names(gs) <- NULL
      return(gs)
    }, loop_hash = loop.obj@loop)
    names(gene_list[[1]]) <- NULL
    shadow_gen_list <- lapply(gene_list[[1]], function(glist)sapply(glist, function(vec)mixedsort(vec)))

    # (genunine) unique gene sets across connectomes
    kpt_agset <- !duplicated(shadow_gen_list)

    # (genuine) remove duplicated gene sets within connectomes
    kpt_wgset <- sapply(gene_list[[1]], function(v_list){
      nv_list <- unique(v_list);
      if (length(nv_list) ==1 ) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })

    stopifnot(length(kpt_agset) == length(kpt_wgset))
    message(sum(!(kpt_agset & kpt_wgset)), " hubs removed due to duplication in idnetical genes contacted within or across hubs.")

    # (genunine) final gene_list
    gene_list[[1]] <- gene_list[[1]][which(kpt_agset & kpt_wgset)]
    ve <- ve[which(kpt_num)][which(kpt_agset & kpt_wgset)]
    ed <- ed[which(kpt_num)][which(kpt_agset & kpt_wgset)]

    # (genuine) gene pairs
    gpair_list <- list()
    gpair_list[[1]] <- rbindlist(lapply(gene_list[[1]], function(glist) {
      gpair <- glist2gpair(glist)$gpair_dat
      return(gpair)
    }))

    if (!random.it) {
      return(gene_list)
    }

    # coordinates of genes
    coord <- rbindlist(lapply(gene_list[[1]], function(gs, info.obj){
        gs <- unique(unlist(gs))
        idx <- chmatch(gs, info.obj@gene[["gene"]])
        stopifnot(all(!is.na(idx)))
        unique(info.obj@gene[idx][, list(chr = chr, start = min(start), end = max(end))])
    }, info.obj = info.obj))

    # (global random) 
    gene_list[[2]] <- coordShulf(coord, info.obj, dout, nshulf = 100, nmin = nmin, nmax = nmax)
    
    # (TAD matched random)
    gene_list[[3]] <- inTADShulf(gene_list[[1]], info.obj)

    for (i in 2:3) {
    	message(ifelse(i==2, "Global random", "TAD matched random"))
    	gpair_list[[i]] <- rbindlist(lapply(gene_list[[i]], function(vec){
        glist2gpair(as.list(vec))$gpair_dat
      }))
    }

    return(gpair_list)
  }
 )