#' @include GRFLoopClass.R

#' @export gene2pairwiseLab
setGeneric(name = "gene2pairwiseLab",
  def = function(gpair_list, info.obj){
    standardGeneric("gene2pairwiseLab")
  }
)

#' @rdname gene2pairwiseLab-methods
setMethod(f = "gene2pairwiseLab",
  signature = c("list", "info"),
  definition = function(gpair_list, info.obj){
    # deg labels for these genes
    col_nm <- colnames(info.obj@gene)[grep("^DEG", colnames(info.obj@gene))]

    deg_list <- lapply(gpair_list, function(dd){
      dd <- as.matrix(dd)
      lapply(1:nrow(dd), function(i) {
        gs <- dd[i,]
        idx <- chmatch(gs, info.obj@gene[["gene"]])
        stopifnot(all(!is.na(idx)))
        deg_dat <- info.obj@gene[idx, col_nm, with = FALSE]
        stopifnot(nrow(deg_dat) == length(gs))
        return(deg_dat)
      })
    })

    lab_list <- lapply(deg_list, function(g_list){
      sapply(g_list, function(dd){
        vec <- dd[[ncol(dd)]]
        if (any(is.na(vec))) {
          return(0)
        } else if (all(vec == "Up")) {
          return(2)
        } else if (all(vec == "Down")) {
          return(1)
        } else if (any(vec == "Up") & any(vec == "Down")) {
          return(-1)
        }
      }) 
    })

    return(lab_list)
  }
)