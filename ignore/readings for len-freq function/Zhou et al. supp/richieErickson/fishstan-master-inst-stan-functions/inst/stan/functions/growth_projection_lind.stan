/* projeciton models for h_growth.stan */

matrix calculate_site_projections_lind_rng(
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
  row_vector[n_sites]  g_growth = // growth coefficient
    exp(mu_beta_cor[2]  + mu_gamma[2]);
  row_vector[n_sites]  t0 = // hypothetical age at which fish's size = 0
    exp(mu_beta_cor[3] + mu_gamma[3]) - 10.0;


  // Use model to make projections
  for(site in 1:n_sites){
    // create projections for each site
    for( p in 1:n_proj){
      site_projections[ site, p] =
        normal_rng(l_inf[site] / (1 + exp(-g_growth[site] *
                   (age_project[p] - t0[site]))),
                   sigma_length);
    }
  }
  return site_projections;
}

vector calculate_hyper_projection_lind_rng(
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
  real g_growth_bar;
  real t0_bar;

  l_inf_bar    = exp(mu_gamma[1]);
  g_growth_bar = exp(mu_gamma[2]);
  t0_bar       = exp(mu_gamma[3]) - 10.0 ;

  // simulate using equations
  for( p in 1:n_proj){
    hyper_projection[p] = normal_rng(
      l_inf_bar / (1 + exp(-g_growth_bar * (age_project[p] - t0_bar))),
      sigma_length);
  }
  return hyper_projection;
}
