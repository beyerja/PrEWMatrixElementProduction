// -----------------------------------------------------------------------------
// Brute forcing this a bit, gonna replace these using bash scripts
static const std::string process = 
  "PROCESS_MARKER"
;

static const std::string output_dir = 
  "OUTPUT_DIR_MARKER"
;

static const int I_n_bins =
  N_BINS_MARKER
;
// -----------------------------------------------------------------------------

void combine_files(
  int energy,
  string I_process_name, 
  double I_xs_LR, 
  double I_xs_RL, 
  double I_xs_LL, 
  double I_xs_RR
) {
  // Convention: I... input variable, O... output variable
  
  //----------------------------------------------------------------------------
  // Input files
  TFile* f_ang = new TFile((output_dir + "/distributions/angular/AngularDistribution_" + to_string(energy) + "GeV_" + I_process_name + ".root").c_str(), "READ");
  TFile* f_tgc = new TFile((output_dir + "/distributions/TGC/AngularTGCDistribution_" + to_string(energy) + "GeV_" + I_process_name + ".root").c_str(), "READ");
  TTree* t_ang = (TTree*)f_ang->Get("AngularDistribution");
  TTree* t_tgc = (TTree*)f_tgc->Get("AngularDistribution");
  
  if ( t_ang->GetEntries() != t_tgc->GetEntries() ) {
    cout << "Oh no, trees don't have same number of entries... Abort!" << endl;
    return;
  }
  
  if ( t_ang->GetEntries() != I_n_bins ) {
    cout << "Number of bins is currently hard coded to 2000 => MUST BE CHANGED IN CODE IF NOT CORRECT!" << endl;
    return;
  }
  
  // Get the values form the trees by connect variables to the tree branches
  const unsigned int I_angle_dim = 3;
  double I_bin_centers[I_angle_dim];
  double I_xs_rel_LR = 0;
  double I_xs_rel_RL = 0;
  double I_xs_rel_RR = 0;
  double I_xs_rel_LL = 0;

  const unsigned int I_n_coefs = 9;
  string * I_tgc_coef_label = new string();		
  double I_tgc_coef_LR[I_n_coefs] {};
  double I_tgc_coef_RL[I_n_coefs] {};
  double I_tgc_coef_RR[I_n_coefs] {};
  double I_tgc_coef_LL[I_n_coefs] {};

  t_ang->SetBranchAddress("angle_center",	I_bin_centers);
  t_ang->SetBranchAddress( "relative_sigma_LR", &I_xs_rel_LR );
  t_ang->SetBranchAddress( "relative_sigma_RL", &I_xs_rel_RL );
  t_ang->SetBranchAddress( "relative_sigma_RR", &I_xs_rel_RR );
  t_ang->SetBranchAddress( "relative_sigma_LL", &I_xs_rel_LL );
  
  t_tgc->SetBranchAddress("PNPC_Label", &I_tgc_coef_label );
  t_tgc->SetBranchAddress("relative_PNPC_LR",	I_tgc_coef_LR	);
  t_tgc->SetBranchAddress("relative_PNPC_RL",	I_tgc_coef_RL	);
  t_tgc->SetBranchAddress("relative_PNPC_RR",	I_tgc_coef_RR	);
  t_tgc->SetBranchAddress("relative_PNPC_LL",	I_tgc_coef_LL	);
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // Counters for normalisation of cross sections
  double sum_I_xs_rel_LR = 0;
  double sum_I_xs_rel_RL = 0;
  double sum_I_xs_rel_RR = 0;
  double sum_I_xs_rel_LL = 0;
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Parameters for RK style tree
  double O_bin_centers[I_n_bins][I_angle_dim] {};

  double O_xs_abs_LR[I_n_bins] {};
  double O_xs_abs_RL[I_n_bins] {};
  double O_xs_abs_RR[I_n_bins] {};
  double O_xs_abs_LL[I_n_bins] {};

  std::string O_tgc_coef_label {};
  double O_tgc_coef_LR[I_n_bins][I_n_coefs] {};
  double O_tgc_coef_RL[I_n_bins][I_n_coefs] {};
  double O_tgc_coef_RR[I_n_bins][I_n_coefs] {};
  double O_tgc_coef_LL[I_n_bins][I_n_coefs] {};
  
  // Combine entries for each bin of distribution into one entry of combined tree
  for (int bin=0; bin<I_n_bins; bin++) {
    t_ang->GetEntry(bin);
    t_tgc->GetEntry(bin);
    
    // Transfer values from input trees to output tree
    for (int a=0; a<I_angle_dim; a++) { O_bin_centers[bin][a] = I_bin_centers[a]; }
    
    O_xs_abs_LR[bin] = I_xs_rel_LR;
    O_xs_abs_RL[bin] = I_xs_rel_RL;
    O_xs_abs_RR[bin] = I_xs_rel_RR;
    O_xs_abs_LL[bin] = I_xs_rel_LL;
    
    if ( bin == 0 ) { O_tgc_coef_label = *I_tgc_coef_label; }
    
    for (int c=0; c<I_n_coefs; c++) { 
      O_tgc_coef_LR[bin][c] = I_tgc_coef_LR[c];
      O_tgc_coef_RL[bin][c] = I_tgc_coef_RL[c];
      O_tgc_coef_RR[bin][c] = I_tgc_coef_RR[c];
      O_tgc_coef_LL[bin][c] = I_tgc_coef_LL[c];
    }
    
    // For normalisation
    sum_I_xs_rel_LR += I_xs_rel_LR;
    sum_I_xs_rel_RL += I_xs_rel_RL;
    sum_I_xs_rel_RR += I_xs_rel_RR;
    sum_I_xs_rel_LL += I_xs_rel_LL;
  } 
  
  // Normalise bins to the total expected cross section
  for (int bin=0; bin<I_n_bins; bin++) {
    if ( sum_I_xs_rel_LR > 0 ) {
      O_xs_abs_LR[bin] *= I_xs_LR / sum_I_xs_rel_LR;
    }
    if ( sum_I_xs_rel_RL > 0 ) {
      O_xs_abs_RL[bin] *= I_xs_RL / sum_I_xs_rel_RL;
    }
    if ( sum_I_xs_rel_RR > 0 ) {
      O_xs_abs_RR[bin] *= I_xs_RR / sum_I_xs_rel_RR;
    }
    if ( sum_I_xs_rel_LL > 0 ) {
      O_xs_abs_LL[bin] *= I_xs_LL / sum_I_xs_rel_LL;
    }
  }
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // New tree in combined file -> same format as RK style files
  TFile* f_com = new TFile((output_dir + "/distributions/combined/Distribution_" + to_string(energy) + "GeV_" + I_process_name + ".root").c_str(), "RECREATE");
  TTree* t_com = new TTree(("MinimizationProcesses" + std::to_string(energy) + "GeV").c_str(), ("MinimizationProcesses" + std::to_string(energy) + "GeV").c_str());
  
  // Values for tree as Robert creates it
  std::string O_process_name = I_process_name;
  int O_n_bins = I_n_bins;
  
  TMatrixT<double> *O_bin_centers_mtx = new TMatrixT<double>(I_n_bins, I_angle_dim);
  TMatrixT<double> *O_tgc_coef_LR_mtx = new TMatrixT<double>(I_n_bins, I_n_coefs);
  TMatrixT<double> *O_tgc_coef_RL_mtx = new TMatrixT<double>(I_n_bins, I_n_coefs);
  TMatrixT<double> *O_tgc_coef_RR_mtx = new TMatrixT<double>(I_n_bins, I_n_coefs);
  TMatrixT<double> *O_tgc_coef_LL_mtx = new TMatrixT<double>(I_n_bins, I_n_coefs);
  O_bin_centers_mtx->Use(I_n_bins, I_angle_dim, &O_bin_centers[0][0]);
  O_tgc_coef_LR_mtx->Use(I_n_bins, I_n_coefs, &O_tgc_coef_LR[0][0]);
  O_tgc_coef_RL_mtx->Use(I_n_bins, I_n_coefs, &O_tgc_coef_RL[0][0]);
  O_tgc_coef_RR_mtx->Use(I_n_bins, I_n_coefs, &O_tgc_coef_RR[0][0]);
  O_tgc_coef_LL_mtx->Use(I_n_bins, I_n_coefs, &O_tgc_coef_LL[0][0]);

  // Branches of RK style tree (Created as Robert created them)
  t_com->Branch("describtion", &O_process_name );
  t_com->Branch("angular_number", &O_n_bins, "angular_number/I");
  t_com->Branch("angular_center", "TMatrixT<double>", &O_bin_centers_mtx, 128000, 0);
  t_com->Branch("differential_sigma_LR", O_xs_abs_LR, "differential_sigma_LR[angular_number]/D" );
  t_com->Branch("differential_sigma_RL", O_xs_abs_RL, "differential_sigma_RL[angular_number]/D" );
  t_com->Branch("differential_sigma_RR", O_xs_abs_RR, "differential_sigma_RR[angular_number]/D" );
  t_com->Branch("differential_sigma_LL", O_xs_abs_LL, "differential_sigma_LL[angular_number]/D" );
  t_com->Branch("differential_PNPC_label", &O_tgc_coef_label );
  t_com->Branch("differential_PNPC_LR", "TMatrixT<double>", &O_tgc_coef_LR_mtx, 256000, 0);
  t_com->Branch("differential_PNPC_RL", "TMatrixT<double>", &O_tgc_coef_RL_mtx, 256000, 0);
  t_com->Branch("differential_PNPC_RR", "TMatrixT<double>", &O_tgc_coef_RR_mtx, 256000, 0);
  t_com->Branch("differential_PNPC_LL", "TMatrixT<double>", &O_tgc_coef_LL_mtx, 256000, 0);
  t_com->Fill();
  
  f_com->Write();
  f_com->Close();
  f_ang->Close();
  f_tgc->Close();
  
  delete O_bin_centers_mtx;
  delete O_tgc_coef_LR_mtx;
  delete O_tgc_coef_RL_mtx;
  delete O_tgc_coef_RR_mtx;
  delete O_tgc_coef_LL_mtx;
}


void CombineDistributions(){
  /** Combine the distributions for the specified process.
  **/
  
  if (process == "ww_sl0muq") {
    combine_files(
      250, // energy,
      "WW_semilep_AntiMuNu", // I_process_name, 
      18781.00 * 0.5 * 0.483, // I_xs_LR, 
      172.73 * 0.5 * 0.483, // I_xs_RL, 
      0, // I_xs_LL, 
      0  // I_xs_RR
    );
    combine_files(
      250, // energy,
      "WW_semilep_MuAntiNu", // I_process_name, 
      18781.00 * 0.5 * 0.483, // I_xs_LR, 
      172.73 * 0.5 * 0.483, // I_xs_RL, 
      0, // I_xs_LL, 
      0  // I_xs_RR
    );
  } else {
    throw std::invalid_argument("Unknown process " + process);
  }
}
