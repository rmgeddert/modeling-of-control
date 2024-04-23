
functions {
  /* Wiener diffusion log-PDF for a single response
  * Args:
  *   y: reaction time data
  *   dec: decision data (0 or 1)
  *   alpha: boundary separation parameter > 0
  *   tau: non-decision time parameter > 0
  *   beta: initial bias parameter in [0, 1]
  *   delta: drift rate parameter
  * Returns:
  *   a scalar to be added to the log posterior
  */
  real wiener_diffusion_lpdf(real y, int dec, real alpha,
                              real tau, real beta, real delta) {
    if (dec == 1) {
      return wiener_lpdf(y | alpha, tau, beta, delta);
    } else {
      return wiener_lpdf(y | alpha, tau, 1 - beta, - delta);
    }
  }
 
  real get_wiener_lpdf_sum(vector RTs, array[] int acc,
                          vector bs, vector ndt, vector bias, vector drift) {
    int N = num_elements(RTs);
    vector[N] lpdf;
    
    for (n in 1:N){
      lpdf[n] = wiener_diffusion_lpdf(RTs[n] | acc[n], bs[n], ndt[n], bias[n], drift[n]);
    }
    
    return sum(lpdf);
  }
}

data {
  int<lower=1> N;  // total number of observations
  vector[N] RTs;  // RTs on each trial
  array[N] int acc;  // accuracy on each trial
  vector<lower=0,upper=1>[N] incongruencies; //if trial is incongruent [1] or congruent [0]
  vector<lower=0,upper=1>[N] switches; //if trial is switch [1] or repeat[1]
  int<lower=1> K; // number of fixed-effects
  vector[N] blockN; // block number (proportion resets when block changes)
  int<lower=1> nSubjects;// total number of subjects
  array[N] int subjects; //which subject does each trial belong to
  int<lower=0, upper=1> priorOnly; //whether to sample from only prior
  int<lower=1> N_pred;  // total number of observations (predicions)
  matrix [N_pred, K] X_pred; // predictions matrix
}

transformed data {
  vector[N] learned_switch_prop = rep_vector(0.0, N);
  vector[N] learned_inc_prop = rep_vector(0.0, N);
  vector[N] averaged_inc_switch_prop = rep_vector(0.0, N);
  matrix[N, K] X; // design matrix
  
  { // block for learning proportions
  print("Learning switch/incongruent proportions...");
  int block_trial_count = 1;
  real cumsum_switch_total = 0;
  real cumsum_inc_total = 0;
  real init_prop_estimate = 0.5;
  
  for (n in 1:N){
    if (n == 1){
      
      // initialize first trial
      learned_switch_prop[n] = init_prop_estimate;
      learned_inc_prop[n] = init_prop_estimate;
      
    } else if (blockN[n] != blockN[n - 1] ) {
      
      // reset when block changes
      learned_switch_prop[n] = init_prop_estimate;
      learned_inc_prop[n] = init_prop_estimate;
      cumsum_switch_total = 0;
      cumsum_inc_total = 0;
      block_trial_count = 1;
      
    } else {
      
      // else increment learned probability
      learned_switch_prop[n] = cumsum_switch_total / block_trial_count;
      learned_inc_prop[n] = cumsum_inc_total / block_trial_count;
      cumsum_switch_total = cumsum_switch_total + switches[n];
      cumsum_inc_total = cumsum_inc_total + incongruencies[n];
      block_trial_count = block_trial_count + 1;
      
    }

  }
  
  // average together switch prop and (reversed) incongruent prop, so that both .75/.75 and .25/.25 average out to .5 switch/inc
  for (n in 1:N){
    averaged_inc_switch_prop[n] = (learned_inc_prop[n] + (1 - learned_switch_prop[n])) / 2;
  }
  
  print("DONE");
  }
  
  { // create design matrix
  print("Creating design matrix...");
  // main effects
  X[, 1] = incongruencies;
  X[, 2] = switches;
  X[, 3] = averaged_inc_switch_prop;
  
  // interactions
  X[, 4] = X[, 1] .* X[, 2]; // congruency vs task sequence
  X[, 5] = X[, 1] .* X[, 3]; // congruency vs averaged_inc_switch_prop
  X[, 6] = X[, 2] .* X[, 3]; // task sequence vs averaged_inc_switch_prop
  X[, 7] = X[, 1] .* X[, 2] .* X[, 3];
  print("DONE");
  }
}

parameters {
  real Intercept_d;  // drift
  real Intercept_bs; // boundary separation
  real Intercept_ndt; // non decision time
  vector[K] b_d;  // drift betas
  vector[K] b_bs;  // boundary separation betas
  vector[K] b_ndt;  // non decsion time betas
  real<lower=0> sd_d;  // drift
  real<lower=0> sd_bs;  // boundary separation
  real<lower=0> sd_ndt;  // non-decision time
  vector[nSubjects] z_d;  // drift subject z scores
  vector[nSubjects] z_bs;  // boundary separation subject z scores
  vector[nSubjects] z_ndt;  // non decision time subject z scores
}

model {
}

generated quantities {
  vector[N_pred] final_drift = rep_vector(0.0, N_pred);
  vector[N_pred] final_bs = rep_vector(0.0, N_pred);
  vector[N_pred] final_ndt = rep_vector(0.0, N_pred);
  vector[N_pred] final_bias = rep_vector(0.5, N_pred);

  // replace zeros with intercept + betas * X_centered
  final_drift += Intercept_d + X_pred * b_d;
  final_bs += Intercept_bs + X_pred * b_bs;
  final_ndt += Intercept_ndt + X_pred * b_ndt;

  // use exp link function for bs and ndt
  final_bs = exp(final_bs);
  final_ndt = exp(final_ndt);
}
