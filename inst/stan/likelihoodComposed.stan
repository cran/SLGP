functions {
  real custom_composed_lpdf(vector epsilon,
  int n,
  int nIntegral,
  int nPredictors,
  int nNeigh,
  matrix functionValues,
  matrix weightMatrix,
  int[,] indMatrix,
  vector weightQuadrature,
  vector trendValues) {
    vector[rows(functionValues)] Z = functionValues * epsilon + trendValues;
    vector[rows(functionValues)] SLGP;

    for (i in 1:nPredictors) {
      int start = 1 + (i - 1) * nIntegral;
      int end = nIntegral + (i - 1) * nIntegral;
      vector[nIntegral] segmentZ = segment(Z, start, nIntegral);

      real maxVal = max(segmentZ);
      segmentZ = exp(segmentZ - maxVal);
      real integralValue = dot_product(segmentZ, weightQuadrature);
      for(j in 1:nIntegral) {
        SLGP[start+j-1] = segmentZ[j]/integralValue;
      }
    }
    vector[n] SLGPcombined;
    int ind;
    for(i in 1:n) {
      SLGPcombined[i] =0;
      for(j in 1:cols(weightMatrix)) {
        ind = indMatrix[i][j];
        SLGPcombined[i] = SLGPcombined[i]  + SLGP[ind]*weightMatrix[i][j];
      }
    }
    return sum(log(SLGPcombined));
  }
}


data {
  int<lower=1> n;
  int<lower=1> nIntegral;
  int<lower=1> nPredictors;
  int<lower=1> nNeigh;
  int<lower=1> p; 										// assuming the length of epsilon
  matrix[nIntegral*nPredictors, p] functionValues;
  matrix[n, nNeigh] weightMatrix;
  int indMatrix[n, nNeigh];
  vector[nIntegral] weightQuadrature;
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
  target += custom_composed_lpdf(epsilon | n, nIntegral, nPredictors, nNeigh, functionValues, weightMatrix, indMatrix, weightQuadrature, trendValues);
}
