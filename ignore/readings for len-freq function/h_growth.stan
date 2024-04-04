// # This originated from https://code.usgs.gov/umesc/quant-ecology/fishstan/-/blob/master/inst/stan/h_growth.stan?ref_type=heads
// # fishStan package 
functions {
#include /functions/growth_models.stan
#include /functions/growth_projection_von_b.stan
#include /functions/growth_projection_von_b_no_t0.stan
#include /functions/growth_projection_gompertz.stan
#include /functions/growth_projection_gq.stan
#include /functions/growth_projection_lind.stan
}
data {
  int<lower = 0> n_fish; // number of fish
  int<lower = 0> n_sites; // number of sites
  int<lower = 0> jj[n_fish] ; // Dummy variable to ID pool
  row_vector<lower = 0>[n_fish] age; // Age of fish
  row_vector<lower = 0>[n_fish] length; // Length of fish
  int<lower = 0, upper = 4> model_type;  // which model 
  // 0 is von_b with t0
  // 1 is von_b without t0
  // 2 is gompertz with t0
  // 3 is gq with t0
  // 4 is lind with t0
  int<lower = 1> n_growth_coef; 
  // 3 for von_b model with t0
  // 2 for von_b model without
  // 3 for bompertz with t0
  // 3 for gq with t0
  // 3 for lind with t0

  // Hyperparameters below:
  real hp_tau;  // Stan's default value is 2.5
  real hp_sigma;  // unsure what is Stan's default value
  real hp_omega; // Stan's default valuedefault 2
  
  // prior for mu_gamma
  real p_mu_gamma;
  real p_mu_gamma_sd;
  int<lower = 0> n_proj; // number of projection points
  vector[n_proj] age_project; // fish ages to simulate lengths for
}
parameters {
  cholesky_factor_corr[n_growth_coef] L_Omega; 
    // prior correlation, Cholesky scale
  matrix[n_growth_coef, n_sites] mu_beta_raw; 
    // This will be transformed into mu_beta
  vector<lower=0>[n_growth_coef] tau; 
    // prior scale
  real<lower = 0> sigma_length; // observation error
  vector[n_growth_coef] mu_gamma; 
    // group coeffs
}
transformed parameters {
  matrix[n_growth_coef, n_sites] mu_beta_cor = 
    make_mu_beta_cor(tau, L_Omega, n_sites, n_growth_coef, mu_beta_raw);

  row_vector[n_fish] growth_model =
      make_growth_model(n_sites, mu_beta_cor, mu_gamma, n_fish, jj,
                        age, model_type);

}
model {
  L_Omega ~ lkj_corr_cholesky(hp_omega);
  to_vector(mu_beta_raw) ~ normal(0, 1);
  tau ~ exponential(1 / hp_tau);
  sigma_length ~ exponential(1 / hp_sigma);
  mu_gamma ~ normal(p_mu_gamma, p_mu_gamma_sd);
  
  length ~ normal(growth_model, sigma_length);
}
generated quantities {

  matrix[n_sites, n_proj] site_projections;
  vector[n_proj] hyper_projection;
  real m_mortality[n_sites];
  real m_mortality_bar;
    
  matrix[n_growth_coef, n_sites] growth_coef;
  row_vector[n_growth_coef] growth_coef_hyper;

  for (ind in 1:n_growth_coef){
    growth_coef[ind] = exp(mu_beta_cor[ind] + mu_gamma[ind]);
    growth_coef_hyper[ind] = exp(mu_gamma[ind]);
  }

  if (model_type == 0) {
    growth_coef[3] += - 10.0;

    m_mortality =
      calculate_m_mortality_von_b(mu_beta_cor, mu_gamma, n_fish, age, n_sites);

    m_mortality_bar = calculate_m_mortality_bar_von_b(mu_beta_cor, mu_gamma);

    site_projections =
      calculate_site_projections_von_b_rng(
              mu_beta_cor, mu_gamma, n_fish, age,
              n_sites, n_proj, age_project, sigma_length);

    hyper_projection =
      calculate_hyper_projection_von_b_rng(
            mu_beta_cor, mu_gamma, n_fish, age, n_sites, n_proj,
            age_project, sigma_length);
  }
  else if (model_type == 1) {
  //  for von b model without t0
    m_mortality =
      calculate_m_mortality_von_b(mu_beta_cor, mu_gamma, n_fish, age, n_sites);

    m_mortality_bar = calculate_m_mortality_bar_von_b(mu_beta_cor, mu_gamma);

    site_projections =
      calculate_site_projections_von_b_no_t0_rng(
              mu_beta_cor, mu_gamma, n_fish,
              age,
              n_sites, n_proj, age_project, sigma_length);


    hyper_projection =
      calculate_hyper_projection_von_b_no_t0_rng(
            mu_beta_cor, mu_gamma, n_fish, age, n_sites, n_proj,
            age_project, sigma_length);
  
  }
  else if (model_type == 2) {
    growth_coef[3] += - 10.0;
  //  for gompertz with t0
    site_projections =
      calculate_site_projections_von_b_no_t0_rng(
              mu_beta_cor, mu_gamma, n_fish, age,
              n_sites, n_proj, age_project, sigma_length);

    hyper_projection =
      calculate_hyper_projection_von_b_no_t0_rng(
            mu_beta_cor, mu_gamma, n_fish, age,
            n_sites, n_proj,
            age_project, sigma_length);
  }
  else if (model_type == 3) {
    growth_coef[3] += - 10.0;
    
    site_projections =
      calculate_site_projections_gq_rng(
              mu_beta_cor, mu_gamma, n_fish, age,
              n_sites, n_proj, age_project, sigma_length);

    hyper_projection =
      calculate_hyper_projection_gq_rng(
            mu_beta_cor, mu_gamma, n_fish, age,
            n_sites, n_proj,
            age_project, sigma_length);
  }
  else if (model_type == 4) {
    growth_coef[3] += - 10.0;
    
    site_projections =
      calculate_site_projections_lind_rng(
              mu_beta_cor, mu_gamma, n_fish, age,
              n_sites, n_proj, age_project, sigma_length);

    hyper_projection =
      calculate_hyper_projection_lind_rng(
            mu_beta_cor, mu_gamma, n_fish, age,
            n_sites, n_proj,
            age_project, sigma_length);
  }
  
}
