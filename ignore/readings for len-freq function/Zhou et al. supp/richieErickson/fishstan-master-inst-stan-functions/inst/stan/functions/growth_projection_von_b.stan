/* projeciton models for h_growth.stan */

real[] calculate_m_mortality_von_b(
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_fish,
  row_vector age,
  int n_sites){
  
    real  m_mortality[n_sites];

    row_vector[n_sites]  l_inf = //theoretical maximum length
      exp(mu_beta_cor[1] + mu_gamma[1]);
    row_vector[n_sites]  k_growth = // growth coefficient
      exp(mu_beta_cor[2]  + mu_gamma[2]);

    // Use model to make projections
    for(site in 1:n_sites){
    // Estimate M from Then et al. 2014
    // * 100 converts from m to cm
    m_mortality[site] =  (4.118 * (k_growth[site] ^ (0.73))) *
        (l_inf[site] * 100) ^(-0.33);
  }

    return m_mortality;
}


matrix calculate_site_projections_von_b_rng(
  matrix mu_beta_cor,
  vector mu_gamma,
  int n_fish,
  row_vector age,
  int n_sites,
  int n_proj,
  vector age_project,
  real sigma_length){

  matrix[n_sites, n_proj] site_projections;

  row_vector[n_sites]  l_inf = //theoretical maximum length
    exp(mu_beta_cor[1] + mu_gamma[1]);
  row_vector[n_sites]  k_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;


  // Use model to make projections
  for(site in 1:n_sites){
    // create projections for each site
    for( p in 1:n_proj){
      site_projections[ site, p] =
        normal_rng(l_inf[site] * (1.0 - exp(-1.0 * k_growth[site] *
                    (age_project[p] - t0[site] ))),
                   sigma_length);
    }
  }
  return site_projections;
}
 
real calculate_m_mortality_bar_von_b(
    matrix mu_beta_cor,
    vector mu_gamma
    ){

  real  m_mortality_bar;

  real l_inf_bar;
  real k_growth_bar;
  real t0_bar;

  l_inf_bar    = exp(mu_gamma[1]);
  k_growth_bar = exp(mu_gamma[2]);

  // Estimate m_mortality using hyperparameters
  m_mortality_bar = 
    (4.118 * (k_growth_bar ^ (0.73))) * (l_inf_bar * 100) ^(-0.33);

  return m_mortality_bar;
}

vector calculate_hyper_projection_von_b_rng(
    matrix mu_beta_cor,
    vector mu_gamma,
    int n_fish,
    row_vector age,
    int n_sites,
    int n_proj,
    vector age_project,
    real sigma_length){

  vector[n_proj] hyper_projection;

  real l_inf_bar;
  real k_growth_bar;
  real t0_bar;

  l_inf_bar    = exp(mu_gamma[1]);
  k_growth_bar = exp(mu_gamma[2]);
  t0_bar       = exp(mu_gamma[3]) - 10.0 ;

  // simulate using equations
  for( p in 1:n_proj){
    hyper_projection[p] = normal_rng(l_inf_bar *
      (1.0 - exp(-1.0 * k_growth_bar * (age_project[p] - t0_bar))),
      sigma_length);
  }
  return hyper_projection;

}
