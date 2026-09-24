functions {
  real custom_composed_lpdf(vector epsilon,
  int n,
  int nIntegral,
  int nPredictors,
  int nNeigh,
  matrix functionValues,
  matrix weightMatrix,
  array[,] int indMatrix,
  vector weightQuadrature,
  vector trendValues) {
    vector[rows(functionValues)] Z = functionValues * epsilon + trendValues;
    vector[rows(functionValues)] logSLGP;

    // Stay on the log scale throughout. The previous version exponentiated,
    // normalised, combined the neighbours on the natural scale and only then
    // took the log: that is both slower (one autodiff node per quadrature
    // node, from the scalar assignment loop) and less stable.
    for (i in 1:nPredictors) {
      int start = 1 + (i - 1) * nIntegral;
      vector[nIntegral] segmentZ = segment(Z, start, nIntegral);
      real maxVal = max(segmentZ);
      real logInt = maxVal
                    + log(dot_product(exp(segmentZ - maxVal), weightQuadrature));
      logSLGP[start:(start + nIntegral - 1)] = segmentZ - logInt;
    }

    vector[n] lp;
    for (i in 1:n) {
      vector[nNeigh] terms;
      for (j in 1:nNeigh)
        terms[j] = log(weightMatrix[i, j] + 1e-300) + logSLGP[indMatrix[i, j]];
      lp[i] = log_sum_exp(terms);
    }
    return sum(lp);
  }
}


data {
  int<lower=1> n;
  int<lower=1> nIntegral;
  int<lower=1> nPredictors;
  int<lower=1> nNeigh;
  int<lower=1> p; 						// assuming the length of epsilon
  matrix[nIntegral*nPredictors, p] functionValues;
  matrix[n, nNeigh] weightMatrix;
  array[n, nNeigh] int indMatrix;
  vector[nIntegral] weightQuadrature;
  real<lower=0> sigma2;						// prior variance of epsilon
  vector[nIntegral*nPredictors] trendValues;
}

transformed data {
  real<lower=0> priorSd = sqrt(sigma2);
}

parameters {
  vector[p] epsilon;
}

model {
  // Prior: always N(0, sigma2 I), so a scalar normal replaces the multivariate
  // one and avoids a p x p Cholesky decomposition at every leapfrog step.
  epsilon ~ normal(0, priorSd);
  // Likelihood:
  target += custom_composed_lpdf(epsilon | n, nIntegral, nPredictors, nNeigh, functionValues, weightMatrix, indMatrix, weightQuadrature, trendValues);
}
