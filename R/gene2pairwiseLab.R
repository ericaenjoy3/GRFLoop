#' @include GRFLoopClass.R

#' @export gene2pairwiseLab
setGeneric(name = "gene2pairwiseLab",
  def = function(gpair_dat, info.obj){
    standardGeneric("gene2pairwiseLab")
  }
)

#' @rdname gene2pairwiseLab-methods
setMethod(f = "gene2pairwiseLab",
  signature = c("data.table", "info"),
  definition = function(gpair_dat, info.obj){

    stopifnot(ncol(gpair_dat) == 2)

    lab_vec <- sapply(transpose(gpair_dat), function(vec){
        glab <- info.obj@gene[chmatch(vec, info.obj@gene[, gene]), DEG_ESC.MEF]
        if (any(is.na(glab))) {
          return(0)
        } else if (all(glab == "Up")) {
          return(2)
        } else if (all(glab == "Down")) {
          return(1)
        } else if (any(glab == "Up") & any(glab == "Down")) {
          return(-1)
        }
    }) 


    return(lab_list)
  }
)