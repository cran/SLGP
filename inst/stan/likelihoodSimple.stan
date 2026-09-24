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
      vector[nIntegral] segmentZ = segment(Z, start, nIntegral);

      real maxVal = max(segmentZ);
      integralValue[i] = maxVal
                         + log(dot_product(exp(segmentZ - maxVal), weightQuadrature));
    }

    return fx1 - dot_product(integralValue, multiplicities);
  }
}


data {
  int<lower=1> n;
  int<lower=1> nIntegral;
  int<lower=1> nPredictors;
  int<lower=1> p; 						// assuming the length of epsilon
  vector[p] meanFvalues;
  matrix[nIntegral*nPredictors, p] functionValues;
  vector[nIntegral] weightQuadrature;
  vector[nPredictors] multiplicities;
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
  target += custom_simple_lpdf(epsilon | n, nIntegral, nPredictors, meanFvalues, functionValues, weightQuadrature, multiplicities, trendValues);
}
