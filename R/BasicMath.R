
## TODO: Do not need examples of everything here

#' Secant Function
#'
#' @param x Input value to take secant
#'
#' @return Secant of x, defined as the reciprocal of the cosine.
#' @export
#'
#' @examples sec(0)
sec <- function(x){
  1/cos(x)
}

#' Expit Function
#'
#' @param x Input value to take expit
#'
#' @return Expit, or the logistic function of x
#' @export
#'
#' @examples expit(0)
expit <- function(x){
  1/(1 + exp(-x))
}

#' Soft Thresholding Function
#'
#' @param x A value which may be a large outlier.
#' @param c Maximum Value which can be achieved.
#'
#' @return Soft threshold of input.  Used to trim outliers.
#' @export
#'
#' @examples SoftThreshold(100)
SoftThreshold <- function(x,c = 10){
  out <- c*tanh(x/c)
  return(out)
}



