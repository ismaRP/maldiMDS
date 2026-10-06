
#' Check whole number
#'
#' @param x `numeric` to check
#' @param tol numreric tolerance. Default is `.Machine$double.eps^0.5`
#' @return FALSE if `x` is not a whole number, TRUE if it is, given a numberic tolerance
is.wholenumber = function(x, tol = .Machine$double.eps^0.5) {
  if (!is.null(x)){
    abs(x - round(x)) < tol
  } else {
    FALSE
  }
}



#' Print only if verbose
#'
#' @param msg Message to print
#' @param verbose `logical` whether to print the message
#'
print_progress = function(msg, verbose) {
  if (verbose) cat(msg)
}



