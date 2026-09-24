#' Predict densities at new covariate locations using a given SLGP model
#'
#' @description
#' `predictSLGP_newNode()` is deprecated; use
#' \code{\link[=predict,SLGP-method]{predict}(object, type = "density")} instead.
#'
#' Computes the posterior predictive probability densities at new covariate
#' points using a fitted Spatial Logistic Gaussian Process (SLGP) model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame containing new covariate values at which to evaluate the SLGP.
#' @param interpolateBasisFun Character string indicating how basis functions are evaluated:
#'   one of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Integer specifying the discretisation step for interpolation (only used if applicable).
#' @param nIntegral Integer specifying the number of quadrature points over the response space.
#' @param normalize Boolean, indicates if we return normalized or unnormalized pdfs. (defaults to TRUE)
#' @param discrete Boolean, indicates if we work with continuous pdfs (default, FALSE) or discrete probabilities
#'
#' @return A data frame combining \code{newNodes} with columns named \code{pdf_1}, \code{pdf_2}, ...,
#' representing the posterior predictive density for each coefficient draw stored in \code{SLGPmodel}.
#'
#' @seealso \code{\link[stats]{predict}} for the recommended user interface.
#'
#' @export
predictSLGP_newNode <- function(SLGPmodel,
                                newNodes,
                                interpolateBasisFun = "WNN",
                                nIntegral=101,
                                nDiscret=101,
                                normalize = TRUE,
                                discrete = FALSE) {
  .Deprecated(msg = paste('predictSLGP_newNode() is deprecated,',
                          'use predict(object, type = "density").'))
  .predict_density(SLGPmodel=SLGPmodel, newNodes=newNodes,
                   interpolateBasisFun = interpolateBasisFun,
                   nIntegral = nIntegral, nDiscret = nDiscret,
                   normalize = normalize, discrete = discrete)
}

## Internal worker: the real body, NO warning, NOT exported.
#' @keywords internal
#' @noRd
.predict_density <- function(SLGPmodel,
                             newNodes,
                             interpolateBasisFun = "WNN",
                             nIntegral = 101,
                             nDiscret = 101,
                             normalize = TRUE,
                             discrete = FALSE) {
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

  if(!discrete){
    quad_w <- rep(1/(nIntegral-1), nIntegral)
    quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
    domain_size <- diff(SLGPmodel@responseRange)
  }else{
    quad_w <- rep(1, nIntegral)
    domain_size <- 1
  }
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
#' @description
#' `predictSLGP_cdf()` is deprecated; use
#' \code{\link[=predict,SLGP-method]{predict}(object, type = "cdf")} instead.
#'
#' Computes posterior predictive cumulative distribution function values at
#' specified response and covariate locations for a fitted Spatial Logistic
#' Gaussian Process (SLGP) model.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame containing both the response column and the
#'   covariate columns at which the CDF is evaluated.
#' @param interpolateBasisFun Character string indicating the interpolation scheme for basis functions:
#'   one of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret discretisation resolution for interpolation (optional).
#' @param nIntegral Number of integration points along the response axis.
#' @param discrete Boolean, indicates if we work with continuous pdfs (default, FALSE) or discrete probabilities
#'
#' @return A data frame combining \code{newNodes} with columns named
#'   \code{cdf_1}, \code{cdf_2}, ..., containing predictive CDF evaluations
#'   for each coefficient draw stored in \code{SLGPmodel}.
#'
#' @seealso \code{\link[stats]{predict}} for the recommended user interface.
#'
#' @export
#'
predictSLGP_cdf <- function(SLGPmodel,
                            newNodes,
                            interpolateBasisFun = "WNN",
                            nIntegral=101,
                            nDiscret=101,
                            discrete=FALSE) {
  .Deprecated(msg = paste('predictSLGP_cdf() is deprecated,',
                          'use predict(object, type = "cdf").'))
  .predict_cdf(SLGPmodel=SLGPmodel, newNodes=newNodes,
               interpolateBasisFun = interpolateBasisFun,
               nIntegral = nIntegral, nDiscret = nDiscret,
               discrete = discrete)
}
## Internal worker: the real body, NO warning, NOT exported.
#' @keywords internal
#' @noRd
.predict_cdf <- function(SLGPmodel,
                         newNodes,
                         interpolateBasisFun = "WNN",
                         nIntegral=101,
                         nDiscret=101,
                         discrete=FALSE) {

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

  if(!discrete){
    quad_w <- rep(1/(nIntegral-1), nIntegral)
    quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
    domain_size <- diff(SLGPmodel@responseRange)
  }else{
    quad_w <- rep(1, nIntegral)
    domain_size <- 1
  }

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
                 intermediateQuantities$weightSamplesToNodesCDF[, j])
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

#' Predict conditional quantiles
#'
#' @description
#' `predictSLGP_quantiles()` is deprecated; use
#' \code{\link[=predict,SLGP-method]{predict}(object, type = "quantiles")} instead.
#'
#' Computes predictive quantiles at specified covariate locations by numerical
#' inversion of the posterior predictive CDF.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame containing the covariate columns at which
#'   quantiles are evaluated.
#' @param probs Numeric vector of probabilities in \eqn{(0, 1)}.
#' @param interpolateBasisFun Character string specifying interpolation scheme:
#'    one of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret discretisation resolution used for CDF inversion.
#' @param nIntegral Number of integration points for computing the SLGP outputs.
#' @param discrete Boolean, indicates if we work with continuous pdfs (default, FALSE) or discrete probabilities
#'
#' @return A data frame containing the covariates in \code{newNodes}, a column
#'   \code{probs}, and columns named \code{qSLGP_1}, \code{qSLGP_2}, ...,
#'   containing predictive quantiles for each coefficient draw stored in
#'   \code{SLGPmodel}.
#'
#' @seealso \code{\link[stats]{predict}} for the recommended user interface.
#'
#' @importFrom stats approx
#'
#' @export
#'
predictSLGP_quantiles <- function(SLGPmodel,
                                  newNodes,
                                  probs,
                                  interpolateBasisFun = "WNN",
                                  nIntegral=101,
                                  nDiscret=101,
                                  discrete=FALSE) {
  .Deprecated(msg = paste('predictSLGP_quantiles() is deprecated',
                          'use predict(object, type = "quantiles").'))
  .predict_quantiles(SLGPmodel=SLGPmodel, newNodes=newNodes,
                     probs=probs, interpolateBasisFun = interpolateBasisFun,
                     nIntegral = nIntegral, nDiscret = nDiscret,
                     discrete = discrete)
}


## Internal worker: the real body, NO warning, NOT exported.
#' @keywords internal
#' @noRd
.predict_quantiles <- function(SLGPmodel,
                               newNodes,
                               probs,
                               interpolateBasisFun = "WNN",
                               nIntegral=101,
                               nDiscret=101,
                               discrete=FALSE) {
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
  if(!discrete){
    quad_w <- rep(1/(nIntegral-1), nIntegral)
    quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
    domain_size <- diff(SLGPmodel@responseRange)
  }else{
    quad_w <- rep(1, nIntegral)
    domain_size <- 1
  }
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


#' Predict conditional moments
#'
#' @description
#' `predictSLGP_moments()` is deprecated; use
#' \code{\link[=predict,SLGP-method]{predict}(object, type = "moments")} instead.
#'
#'
#' Computes raw or centered moments of the posterior predictive distributions
#' at specified covariate locations.
#'
#' @param SLGPmodel An object of class \code{\link{SLGP-class}}.
#' @param newNodes A data frame of new covariate values.
#' @param power Scalar or vector of positive integers indicating the moment orders to compute.
#' @param centered Logical; if \code{TRUE}, computes centered moments. If \code{FALSE}, computes raw moments.
#' @param interpolateBasisFun Interpolation mode for basis functions: \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret discretisation resolution of the response space.
#' @param nIntegral Number of integration points for computing densities.
#' @param discrete Boolean, indicates if we work with continuous pdfs (default, FALSE) or discrete probabilities
#'
#' @return A data frame containing the covariates in \code{newNodes}, a column
#'   \code{power}, and columns named \code{mSLGP_1}, \code{mSLGP_2}, ...,
#'   containing predictive moments for each coefficient draw stored in
#'   \code{SLGPmodel}.
#'
#' @seealso \code{\link[stats]{predict}} for the recommended user interface.
#'
#' @export
#'
predictSLGP_moments <- function(SLGPmodel,
                                newNodes,
                                power,
                                centered=FALSE,
                                interpolateBasisFun = "WNN",
                                nIntegral=101,
                                nDiscret=101,
                                discrete=FALSE) {
  .Deprecated(msg = paste('predictSLGP_moments() is deprecated,',
                          'use predict(object, type = "moments").'))
  .predict_moments(SLGPmodel=SLGPmodel, newNodes=newNodes,
                   power=power, centered=centered,
                   interpolateBasisFun = interpolateBasisFun,
                   nIntegral = nIntegral, nDiscret = nDiscret,
                   discrete = discrete)
}

## Internal worker: the real body, NO warning, NOT exported.
#' @keywords internal
#' @noRd
.predict_moments <- function(SLGPmodel,
                             newNodes,
                             power,
                             centered=FALSE,
                             interpolateBasisFun = "WNN",
                             nIntegral=101,
                             nDiscret=101,
                             discrete=FALSE) {

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
  if(!discrete){
    quad_w <- rep(1/(nIntegral-1), nIntegral)
    quad_w[c(1, nIntegral)] <- quad_w[c(1, nIntegral)]/2
    domain_size <- diff(SLGPmodel@responseRange)
  }else{
    quad_w <- rep(1, nIntegral)
    domain_size <- 1
  }
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

#' Draw conditional samples
#'
#' @description
#' `sampleSLGP()` is deprecated; use
#' \code{\link[=simulate,SLGP-method]{simulate}(object)} instead.
#'
#' Draws samples from the posterior predictive distributions at specified
#' covariate locations using inverse transform sampling based on the estimated
#' predictive CDF.
#'
#' @param SLGPmodel A trained SLGP model object (\code{\link{SLGP-class}}).
#' @param newX A data frame of new covariate values at which to draw samples.
#' @param n Integer or integer vector specifying how many samples to draw at each input point.
#' @param interpolateBasisFun Character string specifying interpolation scheme for basis evaluation.
#'   One of \code{"nothing"}, \code{"NN"}, or \code{"WNN"} (default).
#' @param nDiscret Integer; discretisation step for the response axis.
#' @param nIntegral Integer; number of quadrature points for density approximation.
#' @param seed Optional integer to set a random seed for reproducibility.
#' @param discrete Boolean, indicates if we work with continuous pdfs (default, FALSE) or discrete probabilities
#'
#' @return A data frame containing simulated responses, with the covariate
#'   columns from \code{newX} and one response column named after the response
#'   variable of \code{SLGPmodel}.
#'
#' @seealso \code{\link[stats]{simulate}} for the recommended user interface.
#'
#' @importFrom stats runif approxfun
#'
#' @export
#'
sampleSLGP <- function(SLGPmodel,
                       newX,
                       n,
                       interpolateBasisFun = "WNN",
                       nIntegral=101,
                       nDiscret=101,
                       seed=NULL,
                       discrete=FALSE) {
  .Deprecated(msg = paste('sampleSLGP() is deprecated, use simulate(object).'))
  .simulate_SLGP(SLGPmodel=SLGPmodel, newX=newX, n=n,
                 interpolateBasisFun = interpolateBasisFun,
                 nIntegral=nIntegral, nDiscret=nDiscret, seed=seed,
                 discrete=discrete)
}

## Internal worker: the real body, NO warning, NOT exported.
#' @keywords internal
#' @noRd
.invert_cdf <- function(Fvals, u, r, discrete) {
  if (discrete) {
    u[pmin(findInterval(r, Fvals) + 1L, length(u))]
  } else {
    approx(x = Fvals, y = u, xout = r, ties = "ordered", rule = 2)$y
  }
}

#' @keywords internal
#' @noRd
.simulate_SLGP <- function(SLGPmodel,
                           newX,
                           n,
                           type = c("predictive", "draws"),
                           interpolateBasisFun = "WNN",
                           nIntegral = 101,
                           nDiscret = 101,
                           seed = NULL,
                           discrete = FALSE) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)

  u <- seq(SLGPmodel@responseRange[1], SLGPmodel@responseRange[2], , nIntegral)

  npred <- nrow(newX)
  nsamp <- c(n)
  if (length(nsamp) == 1L) {
    nsamp <- rep(nsamp, npred)
  } else if (length(nsamp) != npred) {
    stop("'n' must have length 1 or nrow(newX).")
  }

  ## Grid at which we want the cdfs
  grid <- expand.grid(u, seq_len(npred))
  grid <- data.frame(cbind(grid[, 1], newX[grid[, 2], ]))
  colnames(grid) <- c(SLGPmodel@responseName, colnames(newX))

  cdfs <- .predict_cdf(SLGPmodel = SLGPmodel,
                       newNodes = grid,
                       interpolateBasisFun = interpolateBasisFun,
                       nIntegral = nIntegral,
                       nDiscret = nDiscret,
                       discrete = discrete)
  cdfs <- as.matrix(cdfs[, -seq_len(ncol(grid)), drop = FALSE])
  ndraws <- ncol(cdfs)

  if (type == "draws" && ndraws == 1L) {
    warning("The model carries a single set of coefficients; ",
            "type = \"draws\" coincides with type = \"predictive\".")
  }
  if (type == "predictive") {
    cdfs <- matrix(rowMeans(cdfs), ncol = 1L)
  }

  res <- lapply(seq_len(npred), function(j) {
    rows <- (j - 1L) * nIntegral + seq_len(nIntegral)
    r <- runif(nsamp[j])
    if (type == "predictive") {
      y <- .invert_cdf(cdfs[rows, 1L], u, r, discrete)
    } else {
      ## one posterior draw per replicate; invert that draw's cdf
      d <- sample.int(ndraws, nsamp[j], replace = TRUE)
      y <- numeric(nsamp[j])
      for (k in unique(d)) {
        sel <- d == k
        y[sel] <- .invert_cdf(cdfs[rows, k], u, r[sel], discrete)
      }
    }
    df <- data.frame(unname(y), unname(newX[rep(j, nsamp[j]), , drop = FALSE]))
    colnames(df) <- c(SLGPmodel@responseName, colnames(newX))
    df
  })
  do.call(rbind, res)
}
