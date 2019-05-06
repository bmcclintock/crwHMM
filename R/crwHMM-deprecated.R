##' Deprecated functions.
##'
##' \code{extract} and \code{pluck} have been deprecated to avoid conflict with raster::extract, tidyr::extract, purrr::pluck. Instead use
##' \code{grab}.
##' @rdname crwHMM-deprecated
##' @param ... ignored
##' @export
##' @rdname crwHMM-deprecated
pluck <- function(...) {
  .Deprecated("grab")
}

##' @export
##' @rdname crwHMM-deprecated
extract <- function(...) {
  .Deprecated("grab")
}
