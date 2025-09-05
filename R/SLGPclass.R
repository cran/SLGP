#' The SLGP S4 Class: Spatial Logistic Gaussian Process Model
#'
#' This S4 class represents a Spatial Logistic Gaussian Process (SLGP) model, designed for
#' modeling conditional or spatially dependent probability distributions. It encapsulates all
#' necessary components for training, sampling, and prediction, including the basis function
#' setup, learned coefficients, and fitted hyperparameters.
#'
#' @aliases SLGP-class
#'
#' @slot formula A \code{formula} specifying the model structure and covariates.
#' @slot data A \code{data.frame} containing the observations used to train the model.
#' @slot responseName A \code{character} string specifying the name of the response variable.
#' @slot covariateName A \code{character} vector specifying the names of the covariates.
#' @slot responseRange A \code{numeric} vector of length 2 indicating the lower and upper bounds of the response.
#' @slot predictorsRange A \code{list} containing:
#'   - \code{predictorsLower}: lower bounds of the covariates;
#'   - \code{predictorsUpper}: upper bounds of the covariates.
#' @slot method A \code{character} string indicating the training method used: one of \{"MCMC", "MAP", "Laplace", "none"\}.
#' @slot p An \code{integer} indicating the number of basis functions used.
#' @slot basisFunctionsUsed A \code{character} string specifying the type of basis functions used:
#'   "inducing points", "RFF", "Discrete FF", "filling FF", or "custom cosines".
#' @slot opts_BasisFun A \code{list} of additional options used to configure the basis functions.
#' @slot BasisFunParam A \code{list} containing the computed parameters of the basis functions,
#'   e.g., Fourier frequencies or interpolation weights.
#' @slot coefficients A \code{matrix} of coefficients for the finite-rank Gaussian process.
#'   Each row corresponds to a realization of the latent field:
#'   \eqn{ Z(x, t) = \sum_{i=1}^p \epsilon_i f_i(x, t) }.
#' @slot hyperparams A \code{list} of hyperparameters, including:
#'   - \code{sigma}: numeric signal standard deviation;
#'   - \code{lengthscale}: a vector of lengthscales for each input dimension.
#' @slot logPost A \code{numeric} value representing the (unnormalized) log-posterior of the model.
#'   Currently available only for MAP and Laplace-trained models.
#'
#' @export
SLGP <- setClass(
  "SLGP",
  slots = c(
    formula = "formula",
    data = "data.frame",
    responseName = "character",
    covariateName = "character",
    responseRange = "numeric",
    predictorsRange = "list",
    method = "character",
    p = "numeric",
    basisFunctionsUsed = "character",
    opts_BasisFun = "list",
    BasisFunParam = "list",
    coefficients = "matrix",
    hyperparams = "list",
    logPost = "numeric"
  )
)
