#' normalize_data: Normalize data to the range \[0, 1\]
#'
#' Scales the response and covariates of a dataset to the unit interval \eqn{[0,1]}.
#' This normalization is required before applying SLGP methods. If range bounds are
#' not provided, they are computed from the data.
#'
#' @param data A data frame containing the dataset.
#' @param predictorNames A character vector of covariate column names.
#' @param responseName A character string specifying the response variable name.
#' @param predictorsUpper Optional numeric vector of upper bounds for covariates.
#' @param predictorsLower Optional numeric vector of lower bounds for covariates.
#' @param responseRange Optional numeric vector of length 2 giving lower and upper bounds for the response.
#'
#' @return A normalized data frame with the same column structure as \code{data}, with values scaled to \eqn{[0,1]}.
#'
#' @keywords internal
#'
normalize_data <- function(data, predictorNames, responseName, predictorsUpper = NULL, predictorsLower = NULL, responseRange = NULL) {
  # Check if predictor and response names exist in the dataframe
  stopifnot("The predictor names are not all in the dataset colnames"=all(predictorNames %in% colnames(data)), responseName %in% colnames(data))

  # Extract predictor and response columns
  predictors <- data[, predictorNames, drop = FALSE]
  response <- data[[responseName]]

  # Use observed values if range information is not provided
  if (is.null(predictorsUpper) || is.null(predictorsLower)) {
    predictorsUpper <- apply(predictors, 2, max, na.rm = TRUE)
    predictorsLower <- apply(predictors, 2, min, na.rm = TRUE)
  }

  if (is.null(responseRange)) {
    responseRange <- c(min(response, na.rm = TRUE), max(response, na.rm = TRUE))
  }

  predictorsScale <- ifelse(abs(predictorsUpper - predictorsLower)<=1e-10, 1, predictorsUpper - predictorsLower)
  # Normalize predictors to the range [0, 1]
  normalized_predictors <- scale(predictors, center = predictorsLower, scale = predictorsScale)

  # Normalize response to the range [0, 1]
  normalized_response <- (response - responseRange[1]) / (responseRange[2] - responseRange[1])

  # Create a new dataframe with normalized values
  normalized_data <- cbind(normalized_response, normalized_predictors)
  colnames(normalized_data) <- c(responseName, predictorNames)

  return(normalized_data)
}

#' pre_comput_nothing: Precompute quantities for SLGP basis evaluation without interpolation
#'
#' Computes intermediate quantities for evaluating basis functions when no interpolation
#' is used. Basis functions are evaluated at the exact covariate and response grid locations.
#'
#' @param normalizedData A data frame with values already normalized to \eqn{[0,1]}.
#' @param predictorNames Character vector of covariate column names.
#' @param responseName Name of the response variable.
#' @param nIntegral Integer, number of points used to discretize the response domain.
#'
#' @return A list of intermediate quantities used in SLGP basis function computation:
#'   \itemize{
#'     \item \code{nodes}: all points where basis functions are evaluated,
#'     \item \code{indNodesToIntegral}: index mapping nodes to response bins,
#'     \item \code{indSamplesToNodes}: index mapping observations to nodes,
#'     \item \code{indSamplesToPredictor}: index mapping observations to unique predictors,
#'     \item \code{weightSamplesToNodes}: interpolation weights (equal to 1 here).
#'   }
#'
#' @keywords internal
#'
pre_comput_nothing <- function(normalizedData, predictorNames, responseName, nIntegral=101) {
  # Extract predictor and response columns
  predictors <- normalizedData[, predictorNames, drop = FALSE]

  # Extract unique predictors and necessary nodes
  uniquePredictors <- unique(predictors)
  nodesIntegral <- expand.grid(seq(0, 1,, nIntegral), seq(nrow(uniquePredictors)))
  nodesIntegral <- cbind(nodesIntegral[, 1], uniquePredictors[nodesIntegral[, 2], ]) #All the values useful in the integral
  colnames(nodesIntegral) <- c(responseName, predictorNames)

  # Where the functions need to be evaluated
  nodes <- unique(rbind(nodesIntegral, normalizedData))

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodesIntegral[, -c(1)])), data.frame(t(uniquePredictors))),
                           rep(NA, nrow(nodes)-nIntegral*nrow(uniquePredictors)))

  # Index of the node to use for each sample point
  indSamplesToNodes <-  as.matrix(match(data.frame(t(normalizedData)), data.frame(t(nodes))))
  # Index of the integral to use for each sample point
  indSamplesToPredictor <-  match(data.frame(t(normalizedData[, -c(1)])), data.frame(t(uniquePredictors)))
  # Corresponding weight
  weightSamplesToNodes <- as.matrix(rep(1, length(indSamplesToNodes)))

  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    indSamplesToPredictor=indSamplesToPredictor,
    weightSamplesToNodes = weightSamplesToNodes
  )

  return(intermediate_quantities)
}

#' pre_comput_NN: Precompute quantities for SLGP basis evaluation with nearest-neighbor interpolation
#'
#' Computes intermediate quantities for evaluating SLGP basis functions using
#' Nearest Neighbor (NN) interpolation over a regular grid in the normalized domain.
#'
#' @param normalizedData A normalized data frame (values in \eqn{[0,1]}).
#' @param predictorNames Character vector of covariate names.
#' @param responseName Name of the response variable.
#' @param nIntegral Number of grid points for discretizing the response domain.
#' @param nDiscret Number of grid points for discretizing the covariate domain.
#'
#' @return A list of intermediate quantities used in SLGP evaluation:
#'   \itemize{
#'     \item \code{nodes}: grid of response × covariates,
#'     \item \code{indNodesToIntegral}: response bin indices,
#'     \item \code{indSamplesToNodes}: sample-to-node index mapping,
#'     \item \code{weightSamplesToNodes}: equal weights for NN interpolation.
#'   }
#'
#' @keywords internal
#'
pre_comput_NN <- function(normalizedData, predictorNames, responseName, nIntegral=101, nDiscret=101) {
  # Lets map the samples to the NN:
  normalizedData[, responseName] <- round(normalizedData[, responseName]*(nIntegral-1))/(nIntegral-1)
  normalizedData[, predictorNames] <- round(normalizedData[, predictorNames]*(nDiscret-1))/(nDiscret-1)

  predictors <- normalizedData[, predictorNames, drop = FALSE]

  # Extract unique predictors and necessary nodes
  uniquePredictors <- unique(predictors)
  # Where the functions need to be evaluated
  nodes <- expand.grid(seq(0, 1,, nIntegral), seq(nrow(uniquePredictors)))
  nodes <- cbind(nodes[, 1], uniquePredictors[nodes[, 2], ]) #All the values useful in the integral
  colnames(nodes) <- c(responseName, predictorNames)

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodes[, -c(1)])), data.frame(t(uniquePredictors))),
                           rep(NA, nrow(nodes)-nIntegral*nrow(uniquePredictors)))

  # Index of the node to use for each sample point
  indSamplesToNodes <-   as.matrix(match(data.frame(t(normalizedData)), data.frame(t(nodes))))
  # Corresponding weight
  weightSamplesToNodes <- as.matrix(rep(1, length(indSamplesToNodes)))

  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    weightSamplesToNodes = weightSamplesToNodes
  )

  return(intermediate_quantities)
}


#' pre_comput_WNN: Precompute quantities for SLGP basis evaluation with weighted nearest-neighbors
#'
#' Computes intermediate quantities for evaluating basis functions via
#' weighted nearest-neighbor (WNN) interpolation on a discretized grid.
#'
#' @param normalizedData Normalized data frame (\eqn{[0,1]}-scaled).
#' @param predictorNames Character vector of covariate names.
#' @param responseName Name of the response variable.
#' @param nIntegral Number of quadrature points for response domain.
#' @param nDiscret Number of discretization steps for covariates.
#'
#' @return A list of intermediate quantities:
#'   \itemize{
#'     \item \code{nodes}: all evaluation points in response × covariates grid,
#'     \item \code{indNodesToIntegral}: indices to map nodes to response bins,
#'     \item \code{indSamplesToNodes}: index mapping from samples to grid nodes,
#'     \item \code{weightSamplesToNodes}: interpolation weights using inverse distance.
#'   }
#'
#' @keywords internal
#'
pre_comput_WNN <- function(normalizedData, predictorNames, responseName,
                           nIntegral=101, nDiscret=101, mode="pdf") {
  # Lets map the samples to the NN:
  normalizedData2<-normalizedData
  normalizedData2[, responseName] <- trunc(normalizedData[, responseName]*(nIntegral-1))/(nIntegral-1)
  normalizedData2[, predictorNames] <- trunc(normalizedData[, predictorNames]*(nDiscret-1))/(nDiscret-1)

  # We need to extract the predictor, and edit it to add all the dimensions
  predictors <- normalizedData2

  neighboursStructure <- t(expand.grid(rep(list(c(0,1/(nDiscret-1))), length(predictorNames)+1)))
  neighboursStructure[1, ] <- neighboursStructure[1, ]*(nDiscret-1)/(nIntegral-1)
  neighboursStructure <- unname(matrix(apply(predictors, 1, FUN=function(x){
    neighboursStructure+x
  }), ncol=ncol(normalizedData2), byrow=TRUE)) # Store all the neighbouring nodes of my points, there may be replicates
  # Or nodes not in the hypercube
  adjacentNodes<- neighboursStructure[apply((neighboursStructure >= 0)*(neighboursStructure <= 1),
                                            1, function(x){prod(x)==1}),  ]
  adjacentNodes<- unique(adjacentNodes)
  # Unique Nodes to be considered
  indSamplesToNodes <- matrix(c(match(data.frame(t(neighboursStructure)), data.frame(t(adjacentNodes)))),
                              ncol=2^(1+length(predictorNames)),
                              byrow=TRUE)
  # Create the weights and adjacency matrix:
  weightSamplesToNodes <- t(sapply(seq(nrow(indSamplesToNodes)), function(i){
    x <-as.numeric(normalizedData[i, ])
    sapply(seq(2^(1+length(predictorNames))), function(j){
      y <- as.numeric(adjacentNodes[indSamplesToNodes[i, j], ])
      return(prod(1-abs(x-y)*c(nIntegral-1, rep(nDiscret-1, length(predictorNames)))))
    })
  }))

  weightSamplesToNodes[is.na(weightSamplesToNodes)] <-0
  indSamplesToNodes[is.na(indSamplesToNodes)]<- 1

  weightSamplesToNodes <- abs(weightSamplesToNodes)
  weightSamplesToNodes <- round(weightSamplesToNodes, 15)
  rm(neighboursStructure, predictors)
  gc()

  # For each of the potential neighbours
  # Extract unique predictors and necessary nodes
  predictors <- round(adjacentNodes[, -c(1), drop=FALSE], 15)
  uniquePredictors <- unique(predictors)
  # Where the functions need to be evaluated
  u <- seq(0, 1,, nIntegral)
  nodes <- expand.grid(u, seq(nrow(uniquePredictors)))
  nodes <- cbind(nodes[, 1], uniquePredictors[nodes[, 2], ]) #All the values useful in the integral
  colnames(nodes) <- c(responseName, predictorNames)

  # Which integral is the node used in
  indNodesToIntegral <-  c(match(data.frame(t(nodes[, -c(1)])), data.frame(t(uniquePredictors))))

  # Index of the node to use for each sample point
  indTemp <- c(match(data.frame(t(adjacentNodes)), data.frame(t(nodes))))
  indSamplesToNodes <- matrix(indTemp[indSamplesToNodes], nrow=nrow(indSamplesToNodes))

  rm(indTemp)


  # Return a list of intermediate quantities
  intermediate_quantities <- list(
    nodes = nodes,
    indNodesToIntegral = indNodesToIntegral,
    indSamplesToNodes = indSamplesToNodes,
    weightSamplesToNodes = weightSamplesToNodes
  )

  if(mode=="cdf"){
    # Need what is the first node of an integral (then all the nodes in an integral are this ind + 1:Nintegral)
    first_idx <- match(unique(indNodesToIntegral), indNodesToIntegral)
    ind_int <- 1:nIntegral

    temp <- sapply(seq(nrow(normalizedData)), function(i){
      indx <- indSamplesToNodes[i, ]
      t <- as.numeric(normalizedData[i, responseName])
      x <- c(normalizedData[i, predictorNames])
      dt <- 1/(nIntegral-1)
      dx <- 1/(nDiscret-1)
      v_weights_t <- ifelse(t<(u-dt), 0,
                            ifelse(t<(u), (t-(u-dt))^2/2/dt,
                                   ifelse(t<(u+dt), dt/2+(t-u)*(u+2*dt-t)/2/dt, dt)))
      v_weights_t[1] <- v_weights_t[1] -dt/2
      v_weights_t <-rep(v_weights_t,length(indx))
      v_weights_x <- rep(0, nIntegral*length(indx))
      v_ind <- rep(ind_int, length(indx))
      uInd <- unique(indNodesToIntegral[indx])
      for(j in seq_along(uInd)){
        y <- c(nodes[first_idx[uInd[j]], predictorNames])
        temp <- abs(x-y)/dx
        v_weights_x[(j-1)*nIntegral+ind_int] <-  prod(ifelse(temp<=1, 1-temp, 0))
        v_ind[(j-1)*nIntegral+ind_int] <- first_idx[uInd[j]]+ind_int-1
      }
      return(c(v_weights_t*v_weights_x, v_ind))
    })
    intermediate_quantities$weightSamplesToNodesCDF <- t(temp[1:(nIntegral*ncol(indSamplesToNodes)), ])
    intermediate_quantities$indSamplesToNodesCDF  <- t(temp[-c(1:(nIntegral*ncol(indSamplesToNodes))), ])
  }


  return(intermediate_quantities)
}
