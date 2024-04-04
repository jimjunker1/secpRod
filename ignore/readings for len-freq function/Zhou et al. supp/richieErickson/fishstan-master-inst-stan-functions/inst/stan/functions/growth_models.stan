/* growth models and functions for h_growth.stan */
  
matrix make_mu_beta_cor(
  vector tau,
  matrix L_Omega,
  int n_sites,
  int n_coef, // 3 for model with t0, 2 for model without
  matrix mu_beta_raw
  ) {
  // correlated site-level variation, without mu_gamma
  matrix[n_coef, n_sites] mu_beta_cor = 
    diag_pre_multiply(tau, L_Omega) * mu_beta_raw;
    // This is part of the MVN non-central parameterization
    // The mean vector (mu_gamma) still needs to be added, but that is done below
  return mu_beta_cor;
  }

row_vector make_von_b_placeholder(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age
  ){
    
  row_vector[n_sites]  l_inf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  k_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
    
  row_vector[n_obs] von_b_placeholder =
    l_inf[jj]  .*  (1.0 - exp( - k_growth[jj] .* (age - t0[jj])));

  return von_b_placeholder;
  }

row_vector make_von_b_no_t0_placeholder(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age
  ){
  row_vector[n_sites]  l_inf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  k_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);

  row_vector[n_obs] von_b_no_t_placeholder =
    l_inf[jj]  .*  (1.0 - exp( - k_growth[jj] .* (age)));

  return von_b_no_t_placeholder;
  }

row_vector make_gompertz_placeholder(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age
  ){
    
  row_vector[n_sites]  l_inf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  g_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
    
  row_vector[n_obs] gompertz_placeholder =
    l_inf[jj]  .*  exp( - exp( - g_growth[jj] .* (age - t0[jj])));

  return gompertz_placeholder;
  }
  
row_vector make_gq_placeholder(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age
  ){
    
  row_vector[n_sites]  w_growth = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  k_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
    
  row_vector[n_obs] gq_placeholder =
    (w_growth[jj] ./ k_growth[jj]) .* 
    (1.0 - exp( - k_growth[jj] .* (age - t0[jj])));

  return gq_placeholder;
  }

row_vector make_lind_placeholder(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age
  ){
  
  row_vector[n_sites]  l_inf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  g_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;
    
  row_vector[n_obs] lind_placeholder =
    l_inf[jj]  ./  (1 + exp(-g_growth[jj] .*(age - t0[jj])) );


  return lind_placeholder;
}

row_vector make_growth_model(
  int n_sites,
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_obs,
  int[] jj,
  row_vector age,
  int model_type){
  
  if (model_type == 0) {
    return make_von_b_placeholder(n_sites, mu_beta_cor, mu_gamma,
                                  n_obs, jj, age);
  } 
  else if (model_type == 1) {
    return make_von_b_no_t0_placeholder(n_sites, mu_beta_cor, mu_gamma,
                                        n_obs, jj, age);
  }
  else if (model_type == 2) {
    return make_gompertz_placeholder(n_sites, mu_beta_cor, mu_gamma,
                                     n_obs, jj, age);
  }
  else if (model_type == 3) {
    return make_gq_placeholder(n_sites, mu_beta_cor, mu_gamma,
                               n_obs, jj, age);
  }
  else if (model_type == 4) {
    return make_lind_placeholder(n_sites, mu_beta_cor, mu_gamma,
                            n_obs, jj, age);
  }
  else {
      // defualt is von_b as well, because Stan requires all if/else to
      // return a value
    return make_von_b_placeholder(n_sites, mu_beta_cor, mu_gamma,
                                  n_obs, jj, age);
    }
  }
