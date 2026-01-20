functions {
  real custom_simple_lpdf(vector epsilon,
  int n,
  int nIntegral,
  int nPredictors,
  vector meanFvalues,
  matrix functionValues,
  vector weightQuadrature,
  vector multiplicities,
  vector trendValues) {
    real fx1 = n * dot_product(epsilon, meanFvalues);

    vector[rows(functionValues)] Z = functionValues * epsilon + trendValues;
    vector[size(multiplicities)] integralValue;

    for (i in 1:size(multiplicities)) {
      int start = 1 + (i - 1) * nIntegral;
      int end = nIntegral + (i - 1) * nIntegral;
      vector[nIntegral] segmentZ = segment(Z, start, nIntegral);

      real maxVal = max(segmentZ);
      segmentZ = exp(segmentZ - maxVal);
      integralValue[i] = maxVal + log(dot_product(segmentZ, weightQuadrature));
    }

    return fx1 - dot_product(integralValue, multiplicities);
  }
}


data {
  int<lower=1> n;
  int<lower=1> nIntegral;
  int<lower=1> nPredictors;
  int<lower=1> p; 										// assuming the length of epsilon
  vector[p] meanFvalues;
  matrix[nIntegral*nPredictors, p] functionValues;
  vector[nIntegral] weightQuadrature;
  vector[nPredictors] multiplicities;
  matrix[p, p] Sigma;      								// Covariance matrix
  vector[p] mean_x;           							// Vector of mean value of the prior
  vector[nIntegral*nPredictors] trendValues;
}

parameters {
  vector[p] epsilon;
}

model {
  // Priors:
  epsilon ~ multi_normal(mean_x, Sigma);
  // Likelihood:
  target += custom_simple_lpdf(epsilon | n, nIntegral, nPredictors, meanFvalues, functionValues, weightQuadrature, multiplicities, trendValues);
}
