
data {
  // data dimensions (N observations, K proteins, P measurements of protein/carbon from Joli et al)
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> P_joli;
  
  // fmol per ug protein measurement means -- mean estimates
  vector[N] cuznsod_meas;
  vector[N] mnsod_meas;
  vector[N] perox2_meas;

  // measurement error of fmol per ug protein measurement -- standard errors
  vector[N] cuznsod_se;
  vector[N] mnsod_se;
  vector[N] perox2_se;

  // proportion of ug protein in Frag -- mean estimates
  vector[N] frag_proportion_meas;
  
  // measurement error of fmol per ug protein measurement -- standard errors
  vector[N] frag_proportion_se;
  
  // covariate data values
  vector[N] temp_vals;
  vector[N] fe_vals;
  vector[N] mn_vals;
  vector[N] light_vals;

  // data on protein per carbon
  vector[P_joli] joli_protein_per_carbon;
  
  // input for prior on protein per carbon;
  real liefer_mean_protein_per_carbon;
  real liefer_sd_protein_per_carbon;
  
  // input for prior on metallation percentage
  int shape1;
  int shape2;
}

parameters {
  // setting up the covariance matrix
  cholesky_factor_corr[K] chol;
  vector<lower=0>[K] sigma;

  // parameters for true observed values
  vector<lower=0>[N] cuznsod_true;
  vector<lower=0>[N] mnsod_true;
  vector<lower=0>[N] perox2_true;

  vector<lower=0, upper=1>[N] frag_proportion_true;// This is a proportion that is therefore restricted b/w 0 and 1
  
  // regression model coefficients
  vector<lower=0>[K] beta_0; // intercept, bounded by zero because all proteins are non-zero values
  vector[K] beta_1_temp;
  vector[K] beta_2_fe;
  vector[K] beta_3_mn;
  vector[K] beta_4_light;
  
  // generated quantities parameters
  real percentage_metallated;
  
  // pars for describing the protein to carbon mass ratio
  real<lower=0> ug_prot_per_ug_carbon_frag;
  real<lower=0> sigma_protein_per_carbon_frag;
}

transformed parameters {
  corr_matrix[K] R; // correlation matrix
  cov_matrix[K] Sigma; // covariance matrix
  matrix[N, K] y_vals; // targeted proteomics values (not normalized by estimated Frag biomass)
  vector[K] y[N]; // collection of values that are estimated fmol protein / ug Fragilariopsis protein. This notation forms an array of K vectors of length N.
  
  // Recomposing the Cholesky factorized matrices to make a covariance matrix
  R = multiply_lower_tri_self_transpose(chol); // R = Lcorr * Lcorr'
  Sigma = quad_form_diag(R, sigma); // quad_form_diag: diag_matrix(sig) * R * diag_matrix(sig)
  
  // Designating columns of the y_vals matrix to be the estimated values for the targeted proteomics data
  y_vals[, 1] = cuznsod_true;
  y_vals[, 2] = mnsod_true;
  y_vals[, 3] = perox2_true;

  // Loop through each element of the matrix and assign it to the array
  for (i in 1:N) {
    for (j in 1:K) {
      // These targeted values are divided by the estimates Frag-derived protein
      y[i, j] = y_vals[i, j]/frag_proportion_true[i]; // divide by the true proportion of Frag protein
    }
  }
}

model {
  //This notation forms an array of K vectors of length N.
  vector[K] mu[N];
  
  frag_proportion_true ~ normal(0.1, 0.5); // weak prior on the proportion of frag in the sample
  
  // measurement error model to estimate the true proportion of Frag protein
  for (i in 1:N){
    frag_proportion_meas[i] ~ normal(frag_proportion_true[i], frag_proportion_se[i]);
  }
  
  // measurement error model for targeted measurements
  for (i in 1:N){
      cuznsod_meas[i] ~ normal(cuznsod_true[i], cuznsod_se[i]);
      mnsod_meas[i] ~ normal(mnsod_true[i], mnsod_se[i]);
      perox2_meas[i] ~ normal(perox2_true[i], perox2_se[i]);
  }
  
  // the means of the multivariate normal are a function of environmental covariates
  for(n in 1:N){
     mu[n] = beta_0 + beta_1_temp*temp_vals[n] + beta_2_fe*fe_vals[n] + beta_3_mn*mn_vals[n] + beta_4_light*light_vals[n];
  }
  
  // y is a matrix of estimated, normalized measurements
  y ~ multi_normal(mu, Sigma);
  
  // priors for gen quantities
  percentage_metallated ~ beta(shape1, shape2);
  
  // prior from Liefer et al data
  ug_prot_per_ug_carbon_frag ~ normal(liefer_mean_protein_per_carbon, liefer_sd_protein_per_carbon);
  
  // estimate the protein per carbon using liefer data as priors
  joli_protein_per_carbon ~ normal(ug_prot_per_ug_carbon_frag, sigma_protein_per_carbon_frag);

}

generated quantities {
  matrix[K,K] Omega;
  matrix[K,K] cov;
  vector[K] x_rand;

  // these are for calculating cofactor per carbon
  real mnfesod_total_metallated;
  real cuznsod_total_metallated;

  real mnfe_per_ug_protein;
  real cuzn_per_ug_protein;

  real mnfe_fmol_per_ug_carbon;
  real cuzn_fmol_per_ug_carbon;

  real mnfe_umol_per_mol_c;
  real cuzn_umol_per_mol_c;

  Omega = multiply_lower_tri_self_transpose(chol);
  cov = quad_form_diag(Omega, sigma);

  // generating the mean value for the multivariate normal
  x_rand = multi_normal_rng(beta_0, cov);

  // transform the MV normal distribution from fmol protein per ug protein

  // first convert to fmol metallated protein per ug protein
  cuznsod_total_metallated = x_rand[1]*percentage_metallated;
  mnfesod_total_metallated = x_rand[2]*percentage_metallated;

  // // convert based on number of cofactors per protein
  mnfe_per_ug_protein = mnfesod_total_metallated*1;
  cuzn_per_ug_protein = cuznsod_total_metallated*1;

  // generate posterior predictions for ug protein per ug carbon;
  real ug_prot_per_ug_carbon_post_pred = normal_rng(ug_prot_per_ug_carbon_frag, sigma_protein_per_carbon_frag);

  // convert to fmol M per C
  mnfe_fmol_per_ug_carbon = mnfe_per_ug_protein*ug_prot_per_ug_carbon_post_pred*1e6;
  cuzn_fmol_per_ug_carbon = cuzn_per_ug_protein*ug_prot_per_ug_carbon_post_pred*1e6;

  // convert to umol M per mol C
  mnfe_umol_per_mol_c = mnfe_fmol_per_ug_carbon*12*1e-9;
  cuzn_umol_per_mol_c = cuzn_fmol_per_ug_carbon*12*1e-9;
}