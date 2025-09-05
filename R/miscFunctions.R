#' Rosenblatt transform to multivariate Student distribution
#'
#' Auxiliary function that maps uniform samples in \[0, 1\]^d to samples from the spectral density
#' of a Matérn kernel (i.e., a multivariate Student distribution).
#'
#' @param x A matrix (or vector) of samples in \[0, 1\]^d to transform.
#' @param dimension Integer. The dimension of the input space.
#' @param MatParam Numeric. The Matérn kernel smoothness parameter (default = 5/2).
#'
#' @return A matrix with transformed coordinates following a multivariate Student distribution.
#'
#' @importFrom stats qt
#'
#' @keywords internal
#'
rosenblatt_transform_multivarStudent <- function(x, dimension, MatParam=5/2){
  ## Starts from a x presumably uniform in [0, 1]. x is either a n*dimension matrix or a dimension-vector
  if(is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if(!is.matrix(x)){
    stop("Please, apply the Rosenblatt transform to a `x` that is either a vector or a matrix")
  }
  if(ncol(x)!=dimension){
    stop("Check `x`'s dimension for the Rosenblatt tranform.")
  }
  # The first marginal is computed:
  newx <- 0*x
  newx[, 1] <- qt(x[, 1], df=2*MatParam)
  partial_sum <- 0
  for(i in seq(2, dimension)){
    #Conditionals
    partial_sum <- partial_sum + x[, i-1]^2
    newx[, i] <- qt(x[, i], df=2*MatParam+i-1)
    newx[, i] <- newx[, i]*sqrt(2*MatParam+partial_sum)/sqrt(2*MatParam +i -1)
  }
  return(newx)
}
