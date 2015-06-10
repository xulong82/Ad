data {
   int<lower=0> N;		# number of outcomes
   int<lower=0> K;		# number of predictors
   matrix<lower=0>[N,K] x;	# predictor matrix
   real y[N];			# outcomes
}

parameters {
   vector<lower=0.001>[K] beta;	# coefficients
   real<lower=0.001>	sigma;	# variation
}

model {
   beta ~ cauchy(0, 10);
   sigma ~ cauchy(0, 10);

   y ~ normal(x * beta, sigma);
}

