#' Predict densities at new covariate locations using a given SLGP model
#'
#' Computes the posterior predictive probability densities at new covariate points
#' using a fitted Spatial Logistic Gaussian Process (SLGP) model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame containing new covariate values at which to evaluate the SLGP.
#' @param interpolateBasisFun Character string indicating how basis functions are evaluated:
#'   one of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Integer specifying the discretization step for interpolation (only used if applicable).
#' @param nIntegral Integer specifying the number of quadrature points over the response space.
#' @param normalise Boolean, indicates if we return normalised or unnormalised pdfs. (defaults to TRUE)
#'
#' @return A data frame combining \code{newNodes} with columns named \code{pdf_1}, \code{pdf_2}, ...,
#' representing the posterior predictive density for each sample of the SLGP.
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#'#' #Create a SLGP model but don't fit it
#' modelPrior <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "none",    # No training
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#'
#' #Let us make 3 draws from the prior
#' nrep <- 3
#' set.seed(8)
#' p <- ncol(modelPrior@coefficients)
#' modelPrior@coefficients <- matrix(rnorm(n=nrep*p), nrow=nrep)
#'
#' # Where to predict the field of pdfs ?
#' dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 5),
#' seq(range_response[1], range_response[2],, 101)))
#' colnames(dfGrid) <- c("age", "medv")
#' predPrior <- predictSLGP_newNode(SLGPmodel=modelPrior,
#'                                  newNodes = dfGrid)
#' }
#' @export
predictSLGP_newNode <- function(SLGPmodel,
                                newNodes,
                                interpolateBasisFun = "WNN",
                                nIntegral=101,
                                nDiscret=101,
                                normalise = TRUE) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  normalizedData <- normalize_data(data=newNodes,
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  trend <- SLGPmodel@trend
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    predictorsUpper <- SLGPmodel@predictorsRange$upper
    predictorsLower <- SLGPmodel@predictorsRange$lower
    responseRange <- SLGPmodel@responseRange
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }

  epsilon <- SLGPmodel@coefficients
  GPvalues <- functionValues %*% t(epsilon) + trendValues
  domain_size <- diff(SLGPmodel@responseRange)

  quad_w <- rep(1/(nIntegral-1), nIntegral)
  quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
  SLGPvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      res<- exp(x-maxval)
      return(res/sum(res*quad_w)/domain_size)
    }))
  })

  res<- sapply(seq(ncol(GPvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  colnames(res) <- paste0("pdf_", seq(ncol(res)))
  res<- cbind(newNodes, res)
  return(res)
}

#' Predict cumulative distribution values at new locations using a SLGP model
#'
#' Computes the posterior cumulative distribution function (CDF) values at specified
#' covariate values using a fitted SLGP model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame with covariate values where the SLGP should be evaluated.
#' @param interpolateBasisFun Character string indicating the interpolation scheme for basis functions:
#'   one of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Discretization resolution for interpolation (optional).
#' @param nIntegral Number of integration points along the response axis.
#'
#' @return A data frame with \code{newNodes} and predicted CDF values, columns named \code{cdf_1}, \code{cdf_2}, ...
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#'#' #Create a SLGP model but don't fit it
#' modelPrior <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "none",    # No training
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#'
#' #Let us make 3 draws from the prior
#' nrep <- 3
#' set.seed(8)
#' p <- ncol(modelPrior@coefficients)
#' modelPrior@coefficients <- matrix(rnorm(n=nrep*p), nrow=nrep)
#'
#' # Where to predict the field of pdfs ?
#' dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 5),
#' seq(range_response[1], range_response[2],, 101)))
#' colnames(dfGrid) <- c("age", "medv")
#' predPriorcdf <- predictSLGP_cdf(SLGPmodel=modelPrior,
#'                                 newNodes = dfGrid)
#' }
#' @export
#'
predictSLGP_cdf <- function(SLGPmodel,
                            newNodes,
                            interpolateBasisFun = "WNN",
                            nIntegral=101,
                            nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  normalizedData <- normalize_data(data=newNodes,
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret,
                                             mode="cdf")
  }
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  trend <- SLGPmodel@trend
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    predictorsUpper <- SLGPmodel@predictorsRange$upper
    predictorsLower <- SLGPmodel@predictorsRange$lower
    responseRange <- SLGPmodel@responseRange
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon) + trendValues

  domain_size <- diff(SLGPmodel@responseRange)
  quad_w <- rep(1/(nIntegral-1), nIntegral)*domain_size
  quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2

  if(interpolateBasisFun =="WNN"){
    SLGPvalues<-sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
        maxval <- max(x)
        res<- exp(x-maxval)
        return(res/sum(res*quad_w))
      }))
    })

    res<- sapply(seq(ncol(SLGPvalues)), function(i){
      unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodesCDF)), function(j){
        return(SLGPvalues[intermediateQuantities$indSamplesToNodesCDF[, j], i]*
                 intermediateQuantities$weightSamplesToNodesCDF[, j]*domain_size)
      }), na.rm = TRUE))
    })

  }else{
    SLGPcvalues<-sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
        maxval <- max(x)
        pdf<- exp(x-maxval)
        cdf<- cumsum(pdf*quad_w)
        cdf<- cdf-min(cdf)
        cdf<- cdf / diff(range(cdf))
        return(cdf)
      }))
    })
    intermediateQuantities$indSamplesToNodes[is.na(intermediateQuantities$indSamplesToNodes)]<- 1
    res<- sapply(seq(ncol(SLGPcvalues)), function(i){
      unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
        return(SLGPcvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
                 intermediateQuantities$weightSamplesToNodes[, j])
      }), na.rm = TRUE))
    })
  }

  colnames(res) <- paste0("cdf_", seq(ncol(res)))
  res<- cbind(newNodes, res)
  return(res)
}

#' Predict quantiles from a SLGP model at new locations
#'
#' Computes quantile values at specified levels (\code{probs}) for new covariate points,
#' based on the posterior CDFs from a trained SLGP model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame of covariate values.
#' @param probs Numeric vector of quantile levels to compute (e.g., 0.1, 0.5, 0.9).
#' @param interpolateBasisFun Character string specifying interpolation scheme: \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Discretization level of the response axis (for CDF inversion).
#' @param nIntegral Number of integration points for computing the SLGP outputs.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item The covariates in \code{newNodes} (repeated per quantile level),
#'     \item A column \code{probs} indicating the quantile level,
#'     \item Columns \code{qSLGP_1}, \code{qSLGP_2}, ... for each posterior sample's quantile estimate.
#'   }
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#' # Train an SLGP model using Laplace estimation and RFF basis
#' modelLaplace <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "Laplace",    # Train using Maximum A Posteriori estimation
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#' dfX <- data.frame(age=seq(range_x[1], range_x[2], 1))
#' # Predict some quantiles, for instance here the first quartile, median, third quartile
#' predQuartiles <- predictSLGP_quantiles(SLGPmodel= modelLaplace,
#'                                        newNodes = dfX,
#'                                        probs=c(0.25, 0.50, 0.75))
#'
#' }
#'
#' @importFrom stats approx
#' @export
predictSLGP_quantiles <- function(SLGPmodel,
                                  newNodes,
                                  probs,
                                  interpolateBasisFun = "WNN",
                                  nIntegral=101,
                                  nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  u <- seq(SLGPmodel@responseRange[1], SLGPmodel@responseRange[2],, nDiscret)
  newNodesX <- newNodes[, predictorNames, drop=FALSE]
  newNodesPred <- expand.grid(u, seq(nrow(newNodesX)))
  IDpred <- newNodesPred[, 2]
  newNodesPred <- as.data.frame(cbind(newNodesPred[, 1],
                                      newNodesX[newNodesPred[, 2], ]))
  colnames(newNodesPred) <- c(responseName, predictorNames)
  normalizedData <- normalize_data(data=newNodesPred[c(responseName, predictorNames)],
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }
  probs <- c(probs)
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  trend <- SLGPmodel@trend
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    predictorsUpper <- SLGPmodel@predictorsRange$upper
    predictorsLower <- SLGPmodel@predictorsRange$lower
    responseRange <- SLGPmodel@responseRange
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }

  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon) + trendValues
  domain_size <- diff(SLGPmodel@responseRange)

  quad_w <- rep(1/(nIntegral-1), nIntegral)
  quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
  SLGPcvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      pdf<- exp(x-maxval)
      cdf<- cumsum(pdf*quad_w)
      cdf<- cdf-min(cdf)
      cdf<- cdf / diff(range(cdf))
      return(cdf)
    }))
  })
  intermediateQuantities$indSamplesToNodes[is.na(intermediateQuantities$indSamplesToNodes)]<- 1
  res<- sapply(seq(ncol(SLGPcvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPcvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  res <- sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(res[, i], IDpred, function(x){
      approx(x=x, y=u, xout=probs, rule=2)$y
    }))
  })
  IDpredX <- c(sapply(seq(nrow(newNodesX)), function(x){rep(x, length(probs))}))
  res<- data.frame(cbind(newNodesX[IDpredX, ], probs, res))
  colnames(res) <- c(predictorNames, "probs", paste0("qSLGP_", seq(ncol(SLGPcvalues))))
  return(res)
}

#' Predict centered or uncentered moments at new locations from a SLGP model
#'
#' Computes statistical moments (e.g., mean, variance, ...) of the posterior predictive
#' distributions at new covariate locations, using a given SLGP model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame of new covariate values.
#' @param power Scalar or vector of positive integers indicating the moment orders to compute.
#' @param centered Logical; if \code{TRUE}, computes centered moments. If \code{FALSE}, computes raw moments.
#' @param interpolateBasisFun Interpolation mode for basis functions: \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Discretization resolution of the response space.
#' @param nIntegral Number of integration points for computing densities.
#'
#' @return A data frame with:
#'   \itemize{
#'     \item Repeated rows of the input covariates,
#'     \item A column \code{power} indicating the moment order,
#'     \item One or more columns \code{mSLGP_1}, \code{mSLGP_2}, ... for the estimated moments across posterior samples.
#'   }
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#' # Train an SLGP model using Laplace estimation and RFF basis
#' modelLaplace <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "Laplace",    # Train using Maximum A Posteriori estimation
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#' dfX <- data.frame(age=seq(range_x[1], range_x[2], 1))
#' predMean <- predictSLGP_moments(SLGPmodel=modelLaplace,
#'                                 newNodes = dfX,
#'                                 power=c(1, 2),
#'                                 centered=FALSE) # Uncentered moments of order 1 and 2
#' predVar <- predictSLGP_moments(SLGPmodel=modelLaplace,
#'                                newNodes = dfX,
#'                                power=c(2),
#'                                centered=TRUE) # Centered moments of order 2 (Variance)
#' }
#' @export
predictSLGP_moments <- function(SLGPmodel,
                                newNodes,
                                power,
                                centered=FALSE,
                                interpolateBasisFun = "WNN",
                                nIntegral=101,
                                nDiscret=101) {
  predictorNames <- SLGPmodel@covariateName
  responseName <-  SLGPmodel@responseName

  u <- seq(SLGPmodel@responseRange[1], SLGPmodel@responseRange[2],, nDiscret)
  newNodesX <- newNodes[, predictorNames, drop=FALSE]
  newNodesPred <- expand.grid(u, seq(nrow(newNodesX)))
  IDpred <- newNodesPred[, 2]
  newNodesPred <- as.data.frame(cbind(newNodesPred[, 1],
                                      newNodesX[newNodesPred[, 2], ]))
  colnames(newNodesPred) <- c(responseName, predictorNames)
  normalizedData <- normalize_data(data=newNodesPred[c(responseName, predictorNames)],
                                   predictorNames=predictorNames,
                                   responseName=responseName,
                                   predictorsUpper = SLGPmodel@predictorsRange$upper,
                                   predictorsLower = SLGPmodel@predictorsRange$lower,
                                   responseRange = SLGPmodel@responseRange)
  # Do we perform exact function evaluation, or we use a grid and interpolate it.
  if(interpolateBasisFun=="nothing"){
    intermediateQuantities <- pre_comput_nothing(normalizedData=normalizedData,
                                                 predictorNames=predictorNames,
                                                 responseName=responseName,
                                                 nIntegral=nIntegral)
  }
  if(interpolateBasisFun =="NN"){
    intermediateQuantities <- pre_comput_NN(normalizedData=normalizedData,
                                            predictorNames=predictorNames,
                                            responseName=responseName,
                                            nIntegral=nIntegral,
                                            nDiscret=nDiscret)
  }
  if(interpolateBasisFun == "WNN"){
    intermediateQuantities <- pre_comput_WNN(normalizedData=normalizedData,
                                             predictorNames=predictorNames,
                                             responseName=responseName,
                                             nIntegral=nIntegral,
                                             nDiscret=nDiscret)
  }
  power <- c(power)
  dimension <- length(predictorNames)+1
  opts_BasisFun <- SLGPmodel@opts_BasisFun
  ## Initialise the basis functions to use
  initBasisFun <- SLGPmodel@BasisFunParam
  lengthscale <- SLGPmodel@hyperparams$lengthscale
  ## Evaluate basis funs on nodes
  functionValues <- evaluate_basis_functions(parameters=initBasisFun,
                                             X=intermediateQuantities$nodes,
                                             lengthscale=lengthscale)
  trend <- SLGPmodel@trend
  if(is.null(trend)){
    trend <- function(df){return(rep(0, nrow(df)))}
    trendValues <- rep(0, nrow(functionValues))
  }else{
    predictorsUpper <- SLGPmodel@predictorsRange$upper
    predictorsLower <- SLGPmodel@predictorsRange$lower
    responseRange <- SLGPmodel@responseRange
    dftrend <- as.data.frame(t(t(as.matrix(intermediateQuantities$nodes))*
                                 c(responseRange[2]-responseRange[1],
                                   predictorsUpper - predictorsLower)+
                                 c(responseRange[1], predictorsLower)))
    colnames(dftrend) <- c(responseName, predictorNames)
    trendValues <- trend(df=dftrend)
    rm(dftrend)
  }
  epsilon <- SLGPmodel@coefficients
  GPvalues <-functionValues %*% t(epsilon)+trendValues
  domain_size <- diff(SLGPmodel@responseRange)

  quad_w <- rep(1/(nIntegral-1), nIntegral)
  quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
  SLGPvalues<-sapply(seq(ncol(GPvalues)), function(i){
    unlist(tapply(GPvalues[, i], intermediateQuantities$indNodesToIntegral, function(x){
      maxval <- max(x)
      res<- exp(x-maxval)
      return(res/sum(res*quad_w)/domain_size)
    }))
  })
  res<- sapply(seq(ncol(GPvalues)), function(i){
    unname(rowSums(sapply(seq(ncol(intermediateQuantities$indSamplesToNodes)), function(j){
      return(SLGPvalues[intermediateQuantities$indSamplesToNodes[, j], i]*
               intermediateQuantities$weightSamplesToNodes[, j])
    }), na.rm = TRUE))
  })
  if(centered){
    res <- sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(res[, i], IDpred, function(x){
        mu <- domain_size*mean(u*x)
        sapply(power, function(y){
          domain_size*mean((u-mu)^y*x)
        })
      }))
    })
  }else{
    res <- sapply(seq(ncol(GPvalues)), function(i){
      unlist(tapply(res[, i], IDpred, function(x){
        sapply(power, function(y){
          domain_size*mean(u^y*x)
        })
      }))
    })
  }

  IDpredX <- c(sapply(seq(nrow(newNodesX)), function(x){rep(x, length(power))}))
  res<- data.frame(cbind(newNodesX[IDpredX, ] ,power, res))
  colnames(res) <- c(predictorNames, "power", paste0("mSLGP_", seq(ncol(SLGPvalues))))
  return(res)
}

#' Draw posterior predictive samples from a SLGP model
#'
#' Samples from the predictive distributions modeled by a SLGP at new covariate inputs.
#' This method uses inverse transform sampling on the estimated posterior CDFs.
#'
#' @param SLGPmodel A trained SLGP model object (\code{\link{SLGP-class}}).
#' @param newX A data frame of new covariate values at which to draw samples.
#' @param n Integer or integer vector specifying how many samples to draw at each input point.
#' @param interpolateBasisFun Character string specifying interpolation scheme for basis evaluation.
#'   One of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Integer; discretization step for the response axis.
#' @param nIntegral Integer; number of quadrature points for density approximation.
#' @param seed Optional integer to set a random seed for reproducibility.
#'
#' @return A data frame containing sampled responses from the SLGP model, with covariate columns from \code{newX}
#' and one response column named after \code{SLGPmodel@responseName}.
#'
#' @examples
#' \donttest{
#' # Load Boston housing dataset
#' library(MASS)
#' data("Boston")
#' # Set input and output ranges manually (you can also use range(Boston$age), etc.)
#' range_x <- c(0, 100)
#' range_response <- c(0, 50)
#'
#' # Train an SLGP model using Laplace estimation and RFF basis
#' modelMAP <- slgp(medv ~ age,        # Use a formula to specify response and covariates
#'                  data = Boston,     # Use the original Boston housing data
#'                  method = "MAP",    # Train using Maximum A Posteriori estimation
#'                  basisFunctionsUsed = "RFF",         # Random Fourier Features
#'                  sigmaEstimationMethod = "heuristic",  # Auto-tune sigma2 (more stable)
#'                  predictorsLower = range_x[1],         # Lower bound for 'age'
#'                  predictorsUpper = range_x[2],         # Upper bound for 'age'
#'                  responseRange = range_response,       # Range for 'medv'
#'                  opts_BasisFun = list(nFreq = 200,     # Use 200 Fourier features
#'                                       MatParam = 5/2), # Matern 5/2 kernel
#'                  seed = 1)                             # Reproducibility
#'
#' # Let's draw new sample points from the SLGP
#'
#' newDataPoints <- sampleSLGP(modelMAP,
#'                             newX = data.frame(age=c(0, 25, 95)),
#'                             n = c(10, 1000, 1), # how many samples to draw at each new x
#'                             interpolateBasisFun = "WNN")
#' }
#' @importFrom stats runif approxfun
#' @export
sampleSLGP <- function(SLGPmodel,
                       newX,
                       n,
                       interpolateBasisFun = "WNN",
                       nIntegral=101,
                       nDiscret=101,
                       seed=NULL) {
  if (!requireNamespace("GoFKernel", quietly = TRUE)) {
    stop("Package 'GoFKernel' could not be used")
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  u <- seq(SLGPmodel@responseRange[1],
           SLGPmodel@responseRange[2],,
           nIntegral)
  # Check if one or many predictors
  npred <- nrow(newX)
  nsamp <- c(n)
  if(length(nsamp)==1 & npred >1){
    nsamp <- rep(nsamp, npred)
  }
  if(length(nsamp)>npred){
    warning("There are more \'n\'s than covariates provided.")
  }
  # Create a grid at which we want the cdfs.
  grid <- expand.grid(u, seq(npred))
  grid <- data.frame(cbind(grid[, 1], newX[grid[, 2], ]))
  colnames(grid)<- c(SLGPmodel@responseName, colnames(newX))


  cdfs <- predictSLGP_cdf(SLGPmodel=SLGPmodel,
                          newNodes=grid,
                          interpolateBasisFun = interpolateBasisFun,
                          nIntegral=nIntegral,
                          nDiscret=nDiscret)
  mean_cdfs <- rowMeans(cdfs[, -c(1:ncol(grid)), drop=FALSE])
  res <- lapply(seq(npred), function(j){
    temp <- mean_cdfs[(j-1)*nIntegral+1:nIntegral]
    f <- approxfun(x=u, y=temp)
    finv<-GoFKernel::inverse(f,
                             lower=SLGPmodel@responseRange[1],
                             upper=SLGPmodel@responseRange[2])
    # plot(f, from=0, to=1)
    # range(temp)
    r <- runif(nsamp[j])
    y<-  sapply(r, finv)
    df <- data.frame(unname(y), unname(newX[rep(j, nsamp[j]), , drop=FALSE]))
    colnames(df)<- c(SLGPmodel@responseName, colnames(newX))
    return(df)
  })
  res<- do.call(rbind, res)
  return(res)
}
