#' @include GRFLoopClass.R

#' @export orderLoop
setGeneric(name = "orderLoop",
  def = function(loop.obj, fet.obj, sm_nm, order_method = c("last_column", "diff", "last_to_first", "quant")){
    standardGeneric("orderLoop")
  }
)

#' @rdname orderLoop-methods
setMethod(f = "orderLoop",
  signature = c("loop", "fet"),
  definition = function(loop.obj, fet.obj, sm_nm, order_method){

    order_method <- match.arg(order_method)
    stopifnot(all(sm_nm %in% fet.obj@hash[["sms"]]))

    order2 <- function(..., decreasing = TRUE){ order(..., decreasing = decreasing) }

    idx <- which(grepl(sm_nm, fet.obj@hash[["sms"]], ignore.case = TRUE) & grepl("Enh", fet.obj@hash[["grps"]], ignore.case = TRUE))

    if (order_method == "last_column") {
      row_idx <- order(fet.obj@dat_list[[idx]][[ncol(fet.obj@dat_list[[idx]])]], decreasing = TRUE)
    }
    if (order_method == "diff") {
      row_idx <- order(fet.obj@dat_list[[idx]][[ncol(fet.obj@dat_list[[idx]])]] - fet.obj@dat_list[[idx]][[1]], decreasing = TRUE)
    }
    if (order_method == "last_to_first") {
      row_idx <- do.call(order2, fet.obj@dat_list[[idx]][,ncol(fet.obj@dat_list[[idx]]):1, with = FALSE])
    }
    if (order_method == "quant") {
      row_idx <- do.call(order2, data.table(quant = as.numeric(loop.obj@split), fet.obj@dat_list[[idx]][,ncol(fet.obj@dat_list[[idx]]):1, with = FALSE]))
    }

    # update split slot of loop object
    loop.obj@split <- factor(as.numeric(loop.obj@split)[row_idx],
      levels = unique(as.numeric(loop.obj@split)[row_idx]))

    # update loop slot of loop object
    loop.obj@loop <- copy(loop.obj@loop[row_idx])
    loop.obj@loop[["rowid"]] <- seq_len(nrow(loop.obj@loop))
    validObject(loop.obj)

    # update dat_list slot of fet object
    fet.obj@dat_list <- lapply(fet.obj@dat_list, function(dat)copy(dat[row_idx]))
    validObject(fet.obj)
    
    return(list(loop.obj = loop.obj, fet.obj = fet.obj))
  }
)
