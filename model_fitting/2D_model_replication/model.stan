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
}

transformed data {
  vector[N] learned_switch_prop = rep_vector(0.0, N);
  vector[N] learned_inc_prop = rep_vector(0.0, N);
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
  print("DONE");
  }
  
  { // create design matrix
  print("Creating design matrix...");
  // main effects
  X[, 1] = incongruencies;
  X[, 2] = switches;
  X[, 3] = learned_switch_prop;
  X[, 4] = learned_inc_prop;
  // 2-way interactions
  X[, 5] = X[, 1] .* X[, 2]; // congruency vs task sequence
  X[, 6] = X[, 1] .* X[, 3]; // congruency vs switch prop
  X[, 7] = X[, 1] .* X[, 4]; // congruency vs inc prop
  X[, 8] = X[, 2] .* X[, 3]; // task sequence vs switch prop
  X[, 9] = X[, 2] .* X[, 4]; // task sequence vs inc prop
  X[, 10] = X[, 3] .* X[, 4]; // switch prop vs inc prop
  X[, 11] = X[, 1] .* X[, 2] .* X[, 3];
  X[, 12] = X[, 1] .* X[, 2] .* X[, 4];
  X[, 13] = X[, 1] .* X[, 3] .* X[, 4];
  X[, 14] = X[, 2] .* X[, 3] .* X[, 4];
  X[, 15] = X[, 1] .* X[, 2] .* X[, 3] .* X[, 4];
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
  if (!priorOnly){
    
    // initialize linear predictor terms as zeros
    vector[N] drift = rep_vector(0.0, N);
    vector[N] bs = rep_vector(0.0, N);
    vector[N] ndt = rep_vector(0.0, N);
    vector[N] bias = rep_vector(0.5, N);
    
    // replace zeros with intercept
    drift += Intercept_d + X * b_d;
    bs += Intercept_bs + X * b_bs;
    ndt += Intercept_ndt + X * b_ndt;
    
    // add random effects
    for (n in 1:N) {
      drift[n] += sd_d * z_d[subjects[n]];
      bs[n] += sd_bs * z_bs[subjects[n]];
      ndt[n] += sd_ndt * z_ndt[subjects[n]];
    }
    
    // convert from log scale
    bs = exp(bs);
    ndt = exp(ndt);
  
    target += get_wiener_lpdf_sum(RTs, acc, bs, ndt, bias, drift);
  }
  
  // priors
  target += normal_lpdf(Intercept_d | 2, 0.5);
  target += normal_lpdf(Intercept_bs | .5, .2);
  target += normal_lpdf(Intercept_ndt | -2, 0.6);
  target += normal_lpdf(b_d | 0, 0.2);
  target += normal_lpdf(b_bs | 0, 0.2);
  target += normal_lpdf(b_ndt | 0, 0.2);
  target += student_t_lpdf(sd_d | 3, 0, 1)
    - 1 * student_t_lccdf(0 | 3, 0, 1);
  target += student_t_lpdf(sd_ndt | 3, 0, 0.25)
     - 1 * student_t_lccdf(0 | 3, 0, 0.25);
  target += student_t_lpdf(sd_bs | 3, 0, 0.25)
    - 1 * student_t_lccdf(0 | 3, 0, 0.25);
  target += std_normal_lpdf(z_d);
  target += std_normal_lpdf(z_ndt);
  target += std_normal_lpdf(z_bs);
}

generated quantities {
  //save out log likelihood
  vector[N] log_lik;

  if (!priorOnly){

    vector[N] final_drift = rep_vector(0.0, N);
    vector[N] final_bs = rep_vector(0.0, N);
    vector[N] final_ndt = rep_vector(0.0, N);
    vector[N] final_bias = rep_vector(0.5, N);

    // replace zeros with intercept + betas * X_centered
    final_drift += Intercept_d + X * b_d;
    final_bs += Intercept_bs + X * b_bs;
    final_ndt += Intercept_ndt + X * b_ndt;

    // add in random effects terms for drift, bs, and ndt
    for (n in 1:N) {
      final_drift[n] += sd_d * z_d[subjects[n]];
      final_bs[n] += sd_bs * z_bs[subjects[n]];
      final_ndt[n] += sd_ndt * z_ndt[subjects[n]];
    }

    // use exp link function for bs and ndt
    final_bs = exp(final_bs);
    final_ndt = exp(final_ndt);

    for (n in 1:N){
      log_lik[n] =  wiener_diffusion_lpdf(RTs[n] | acc[n], final_bs[n], final_ndt[n], final_bias[n], final_drift[n]);
    }

  }
}
