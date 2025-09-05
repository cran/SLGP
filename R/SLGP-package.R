#' SLGP: A package for spatially dependent probability distributions
#'
#' The `SLGP` package implements Spatial Logistic Gaussian Processes (SLGP) for the flexible modeling
#' of conditional and spatially dependent probability distributions. The SLGP framework leverages
#' basis-function expansions and sample-based inference (e.g., MAP, Laplace, MCMC) for efficient
#' density estimation and uncertainty quantification. This package includes functionality to define,
#' train, and sample from SLGP models, as well as visualization and diagnostic tools.
#'
#' @section SLGP functions:
#' The core functions in the package include:
#' - \code{\link{slgp}}: trains an SLGP model from formula, data, and hyperparameters.
#' - \code{\link{predictSLGP_moments}}: computes posterior predictive means and variances.
#' - \code{\link{predictSLGP_quantiles}}: computes posterior predictive quantiles.
#' - \code{\link{sampleSLGP}}: draws samples from the posterior predictive SLGP.
#' - \code{\link{retrainSLGP}}: retrains a fitted SLGP object with new parameters or method.
#'
#' @name SLGP-package
#' @docType package
#' @useDynLib SLGP, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Gautier, Athénaïs (2023). "Modelling and Predicting Distribution-Valued Fields with Applications to Inversion Under Uncertainty." Thesis, Universität Bern, Bern.
#' See the thesis online at \url{https://boristheses.unibe.ch/4377/}
#'
#' @keywords internal
"_PACKAGE"
