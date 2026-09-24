## ---------------------------------------------------------------------------
## Analytic Hessian of the negative log-posterior
##
## rstan::optimizing(hessian = TRUE) differentiates the gradient by finite
## differences, which costs O(p) extra gradient evaluations *after* the
## optimiser has converged. For p in the hundreds this dominates the whole fit.
##
## The Hessian is available in closed form (see the appendix of the companion
## paper): writing Y_k for a random variable with the fitted SLGP density at the
## k-th distinct covariate value,
##
##     d2 l / de_i de_j  =  sum_k n_k Cov( f_i(x_k, Y_k), f_j(x_k, Y_k) ),
##
## to which the N(0, sigma2 I) prior contributes diag(1 / sigma2). The prior
## term is also what makes the result strictly positive definite, hence the
## uniqueness of the mode.
## ---------------------------------------------------------------------------

## Per-block mean and covariance of the basis functions under the fitted
## density. Returns the list of covariance matrices only if 'want_cov' is TRUE,
## since the WNN branch needs them one at a time.
#' @keywords internal
#' @noRd
.slgp_block_moments <- function(functionValues, epsilon, trendValues,
                                nIntegral, nBlocks, weightQuadrature) {
  z <- as.vector(functionValues %*% epsilon) + trendValues
  dens <- numeric(length(z))
  means <- matrix(0, nBlocks, ncol(functionValues))
  for (k in seq_len(nBlocks)) {
    idx <- (k - 1L) * nIntegral + seq_len(nIntegral)
    zk <- z[idx]
    ek <- exp(zk - max(zk))
    ek <- ek / sum(ek * weightQuadrature)          # density at the nodes
    dens[idx] <- ek
    means[k, ] <- crossprod(functionValues[idx, , drop = FALSE],
                            ek * weightQuadrature)
  }
  list(dens = dens, means = means)
}

## Covariance of the basis functions under the density of one block.
#' @keywords internal
#' @noRd
.slgp_block_cov <- function(functionValues, idx, dens, mean_k, weightQuadrature) {
  Fk <- functionValues[idx, , drop = FALSE]
  w  <- dens[idx] * weightQuadrature
  crossprod(Fk, Fk * w) - tcrossprod(mean_k)
}

## Hessian for the "nothing" / "NN" likelihoods, where the linear term in
## epsilon contributes nothing and only the K normalising integrals matter.
#' @keywords internal
#' @noRd
.slgp_hessian_simple <- function(functionValues, epsilon, trendValues,
                                 nIntegral, multiplicities, weightQuadrature,
                                 sigma2) {
  nBlocks <- length(multiplicities)
  mom <- .slgp_block_moments(functionValues, epsilon, trendValues,
                             nIntegral, nBlocks, weightQuadrature)
  H <- matrix(0, ncol(functionValues), ncol(functionValues))
  for (k in seq_len(nBlocks)) {
    if (multiplicities[k] == 0) next
    idx <- (k - 1L) * nIntegral + seq_len(nIntegral)
    H <- H + multiplicities[k] *
      .slgp_block_cov(functionValues, idx, mom$dens, mom$means[k, ],
                      weightQuadrature)
  }
  H + diag(1 / sigma2, ncol(functionValues))
}

## Hessian for the "WNN" likelihood, where each observation contributes
## log( sum_j w_ij p(t_ij | x_ij) ) and the neighbours do not separate.
#' @keywords internal
#' @noRd
.slgp_hessian_wnn <- function(functionValues, epsilon, trendValues,
                              nIntegral, nBlocks, weightMatrix, indMatrix,
                              weightQuadrature, sigma2) {
  p <- ncol(functionValues)
  mom <- .slgp_block_moments(functionValues, epsilon, trendValues,
                             nIntegral, nBlocks, weightQuadrature)
  n <- nrow(indMatrix)

  ## a_ij = w_ij * density at the neighbour node; L_i = sum_j a_ij
  a <- weightMatrix * matrix(mom$dens[indMatrix], nrow = n)
  Li <- rowSums(a)
  Li[Li <= 0] <- .Machine$double.xmin
  ab <- a / Li                                   # a_ij / L_i

  blockOf <- ((indMatrix - 1L) %/% nIntegral) + 1L

  ## Term 1: - sum_i sum_j ab_ij ( v_ij v_ij' ), with v_ij = F_ij - m_block
  V <- functionValues[as.vector(indMatrix), , drop = FALSE] -
       mom$means[as.vector(blockOf), , drop = FALSE]
  wv <- as.vector(ab)
  H <- -crossprod(V, V * wv)

  ## Term 2: + sum_i g_i g_i', with g_i = sum_j ab_ij v_ij
  G <- matrix(0, n, p)
  for (j in seq_len(ncol(indMatrix))) {
    sl <- (j - 1L) * n + seq_len(n)
    G <- G + V[sl, , drop = FALSE] * ab[, j]
  }
  H <- H + crossprod(G)

  ## Term 3: + sum_b s_b C_b, grouping the block covariances
  s <- tapply(as.vector(ab), as.vector(blockOf), sum)
  for (b in names(s)) {
    k <- as.integer(b)
    idx <- (k - 1L) * nIntegral + seq_len(nIntegral)
    H <- H + s[[b]] * .slgp_block_cov(functionValues, idx, mom$dens,
                                      mom$means[k, ], weightQuadrature)
  }
  H + diag(1 / sigma2, p)
}

## Dispatcher: builds the Hessian of the negative log-posterior at 'epsilon'
## from the same data list that was handed to Stan.
#' @keywords internal
#' @noRd
.slgp_hessian <- function(stan_data, epsilon, interpolateBasisFun, sigma2) {
  fv <- stan_data$functionValues
  tv <- stan_data$trendValues
  if (length(tv) != nrow(fv)) tv <- rep(0, nrow(fv))
  if (identical(interpolateBasisFun, "WNN")) {
    .slgp_hessian_wnn(fv, epsilon, tv,
                      nIntegral = stan_data$nIntegral,
                      nBlocks = stan_data$nPredictors,
                      weightMatrix = stan_data$weightMatrix,
                      indMatrix = stan_data$indMatrix,
                      weightQuadrature = stan_data$weightQuadrature,
                      sigma2 = sigma2)
  } else {
    .slgp_hessian_simple(fv, epsilon, tv,
                         nIntegral = stan_data$nIntegral,
                         multiplicities = stan_data$multiplicities,
                         weightQuadrature = stan_data$weightQuadrature,
                         sigma2 = sigma2)
  }
}
