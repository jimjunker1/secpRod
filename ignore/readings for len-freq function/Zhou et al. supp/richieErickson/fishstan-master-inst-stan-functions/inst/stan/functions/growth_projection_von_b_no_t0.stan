/* projeciton models for h_growth.stan
 * most functiosn are in grwoth_projectio_von_b.stan.
 * This file only contains functiosn that have had t0 removed.
 */
matrix calculate_site_projections_von_b_no_t0_rng(
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



  // Use model to make projections
  for(site in 1:n_sites){
    // create projections for each site
    for( p in 1:n_proj){
      site_projections[ site, p] =
        normal_rng(l_inf[site] * (1.0 - exp(-1.0 * k_growth[site] *
                    (age_project[p] ))),
                   sigma_length);
    }
  }
  return site_projections;
}
 
vector calculate_hyper_projection_von_b_no_t0_rng(
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

  l_inf_bar    = exp(mu_gamma[1]);
  k_growth_bar = exp(mu_gamma[2]);

  // simulate using equations
  for (p in 1:n_proj) {
    hyper_projection[p] = normal_rng(l_inf_bar *
      (1.0 - exp(-1.0 * k_growth_bar * (age_project[p]))),
      sigma_length);
  }
  return hyper_projection;

}
