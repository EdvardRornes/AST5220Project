#include "../include/Perturbations.h"

//=============
// Constructors
//=============

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//===================
// Do all the solving
//===================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//=====================================
// The main work: integrate all the 
// perturbations and spline the results
//=====================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");
  // Set up x- and k-arrays with the latter having log spacing
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array     = exp(log_k_array);
  Vector x_array     = Utils::linspace(x_start, x_end, n_x);

  // Arrays to fill in with relevant solution later
  Vector vec_delta_cdm(n_k*n_x,0);
  Vector vec_delta_b(n_k*n_x,0);
  Vector vec_v_cdm(n_k*n_x,0);
  Vector vec_v_b(n_k*n_x,0);
  Vector vec_Psi(n_k*n_x,0);
  Vector vec_Phi(n_k*n_x,0);
  std::vector<Vector> vec_Theta(Constants.n_ell_theta, Vector(n_k*n_x,0));

  // Loop over all wavenumbers
  #pragma omp parallel for schedule(dynamic, 1) // Parallel computing
  for(int ik = 0; ik < n_k; ik++){
    // Progress bar
    if ((10*ik)/n_k != (10*ik + 10)/n_k) {
      std::cout << (100*ik + 100)/n_k << "% " << std::flush;
      if (ik == n_k - 1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    int idx_end_tc  = get_tight_coupling_time_idx(k, x_array);
    double x_end_tc = x_array[idx_end_tc];

    int points_in_tc   = idx_end_tc + 1;
    int points_in_full = n_x - points_in_tc;

    // Make array for tight-coupling
    Vector x_tc_array(points_in_tc);
    for (int i = 0; i < points_in_tc; i++){
      x_tc_array[i] = x_array[i];
    }

    // Make array for full system
    Vector x_full_array(points_in_full);
    for (int i = 0; i < points_in_full; i++){
      x_full_array[i] = x_array[points_in_tc + i];
    }

    // The tight coupling ODE system
    ODEFunction dydx_tc = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };
    ODESolver ODE_tc;
    // Set up initial conditions in the tight coupling regime
    auto y_tc_ini = set_ic(x_start, k);
    ODE_tc.solve(dydx_tc, x_tc_array, y_tc_ini);

    auto sol_tc = ODE_tc.get_data();
  

    // Fill in tight coupling data
    for (int ix = 0; ix < points_in_tc; ix++){
      int idx = ix + n_x*ik;
      auto y_tc = sol_tc[ix];

      // tc quantities
      double &delta_cdm =  y_tc[Constants.ind_deltacdm_tc];
      double &delta_b   =  y_tc[Constants.ind_deltab_tc];
      double &v_cdm     =  y_tc[Constants.ind_vcdm_tc];
      double &v_b       =  y_tc[Constants.ind_vb_tc];
      double &Phi       =  y_tc[Constants.ind_Phi_tc];
      double *Theta     = &y_tc[Constants.ind_start_theta_tc];

      double x = x_tc_array[ix];

      // Constants, cosmo and recombination values
      const int n_ell_theta = Constants.n_ell_theta;
      const double c = Constants.c;
      double Hp      = cosmo->Hp_of_x(x);
      double H0      = cosmo->get_H0();
      double OmegaR0 = cosmo->get_OmegaR(0);
      double dtaudx  = rec->dtaudx_of_x(x);
      double ck_Hp   = c*k/Hp;

      vec_delta_cdm[idx] = delta_cdm;
      vec_delta_b[idx]   = delta_b;
      vec_v_cdm[idx]     = v_cdm;
      vec_v_b[idx]       = v_b;
      vec_Phi[idx]       = Phi;
      vec_Psi[idx]       = -Phi - 12.0*pow(H0/(c*k*exp(x)), 2.0)*OmegaR0*(-20.0*ck_Hp/(45.0*dtaudx)*Theta[1]);

      vec_Theta[0][idx] = Theta[0];
      vec_Theta[1][idx] = Theta[1];
      vec_Theta[2][idx] = -20.0/45.0*ck_Hp/dtaudx*Theta[1];
      for (int ell = 3; ell < n_ell_theta; ell++){
        vec_Theta[ell][idx] = -ell/(2.0*ell + 1)*ck_Hp/dtaudx*vec_Theta[ell - 1][idx];
      } 
    }
    
    
    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    ODESolver ODE_full;
    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    Vector y_full_ini = sol_tc[idx_end_tc];
    auto y_full_ini_vec = set_ic_after_tight_coupling(y_full_ini, x_end_tc, k);
    ODE_full.solve(dydx_full, x_full_array, y_full_ini_vec);
    auto sol_full = ODE_full.get_data();

    // Fill in full data
    for (int ix = points_in_tc; ix < n_x; ix++){
      int idx = ix + n_x*ik;
      int real_ix = ix - points_in_tc;
      auto y = sol_full[real_ix];

      // Full quantities
      double &delta_cdm =  y[Constants.ind_deltacdm];
      double &delta_b   =  y[Constants.ind_deltab];
      double &v_cdm     =  y[Constants.ind_vcdm];
      double &v_b       =  y[Constants.ind_vb];
      double &Phi       =  y[Constants.ind_Phi];
      double *Theta     = &y[Constants.ind_start_theta];

      double x = x_full_array[real_ix];

      // Constants, cosmo and recombination values
      const int n_ell_theta = Constants.n_ell_theta;
      const double c = Constants.c;
      double H0      = cosmo->get_H0();
      double OmegaR0 = cosmo->get_OmegaR(0);
      double dtaudx  = rec->dtaudx_of_x(x);

      vec_delta_cdm[idx] = delta_cdm;
      vec_delta_b[idx]   = delta_b;
      vec_v_cdm[idx]     = v_cdm;
      vec_v_b[idx]       = v_b;
      vec_Phi[idx]       = Phi;
      vec_Psi[idx]       = -Phi - 12.0*pow(H0/(c*k*exp(x)), 2.0)*OmegaR0*Theta[2];

      for (int ell = 0; ell < n_ell_theta; ell++){
        vec_Theta[ell][idx] = Theta[ell];
      }
    }
    
  }
  Utils::EndTiming("integrateperturbation");

  //===============
  // Spline results
  //===============
  Phi_spline.create(x_array, k_array, vec_Phi);
  Psi_spline.create(x_array, k_array, vec_Psi);
  delta_cdm_spline.create(x_array, k_array, vec_delta_cdm);
  delta_b_spline.create(x_array, k_array, vec_delta_b);
  v_cdm_spline.create(x_array, k_array, vec_v_cdm);
  v_b_spline.create(x_array, k_array, vec_v_b);
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int ell = 0; ell < Constants.n_ell_theta; ell++){
    Theta_spline[ell].create(x_array, k_array, vec_Theta[ell]);
  }
}

//========================================================
// Set the initial conditions in the tight coupling regime
//========================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // References to the tight coupling quantities
  double &delta_cdm =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b   =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm     =  y_tc[Constants.ind_vcdm_tc];
  double &v_b       =  y_tc[Constants.ind_vb_tc];
  double &Phi       =  y_tc[Constants.ind_Phi_tc];
  double *Theta     = &y_tc[Constants.ind_start_theta_tc];

  // Constants, cosmo and recombination values
  const double c = Constants.c;
  double Hp = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);
  double ck_Hp = c*k/Hp;
  double Psi = -2.0/3.0;

  // Scalar quantities
  Phi        = -Psi;
  delta_cdm  = -3.0/2.0*Psi;
  delta_b    = -3.0/2.0*Psi;
  v_cdm      = -1.0/2.0*ck_Hp*Psi;
  v_b        = -1.0/2.0*ck_Hp*Psi;

  // Photon temperature perturbations (Theta_ell)
  Theta[0] = -Psi/2.0;
  Theta[1] = 1.0/6.0*ck_Hp*Psi;

  return y_tc;
}

//================================================================
// Set IC for the full ODE system after tight coupling regime ends
//================================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================
  // Compute where in the y array each component belongs and 
  // where corresponding components are located in the y_tc array
  //=============================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta      = Constants.n_ell_theta;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc   = Constants.n_ell_theta_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc   =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc     =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc       =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc       =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc     = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm =  y[Constants.ind_deltacdm_tc];
  double &delta_b   =  y[Constants.ind_deltab_tc];
  double &v_cdm     =  y[Constants.ind_vcdm_tc];
  double &v_b       =  y[Constants.ind_vb_tc];
  double &Phi       =  y[Constants.ind_Phi_tc];
  double *Theta     = &y[Constants.ind_start_theta_tc];

  //==================================================================
  // Fill in the initial conditions for the full equation system below
  //==================================================================

  // Constants, cosmo and recombination values
  const double c = Constants.c;
  double Hp      = cosmo->Hp_of_x(x);
  double dtaudx  = rec->dtaudx_of_x(x);
  double ck_Hp   = c*k/Hp;

  // Scalar quantities (Gravitational potental, baryons and CDM)
  Phi       = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;

  // Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = Theta_tc[2];
  for (int ell = 3; ell < n_ell_theta; ell++){
    Theta[ell] = -(ell + 0.0)/(2.0*ell + 1.0)*ck_Hp/dtaudx*Theta_tc[ell-1];
  }

  return y;
}

//=================================
// The time when tight coupling end
//=================================

double Perturbations::get_tight_coupling_time_idx(const double k, Vector x_arr) const{
  //==================================================
  // Compute and return x for when tight coupling ends
  //==================================================
  int idx_tc = 0;
  bool tc = true;
  double x_current;
  double ck_Hp;
  double dtaudx;
  while (tc && idx_tc < x_arr.size()) {
    x_current = x_arr[idx_tc];
    ck_Hp = Constants.c*k / cosmo->Hp_of_x(x_current);
    dtaudx = abs(rec->dtaudx_of_x(x_current));
    if (abs(dtaudx) < 10.0 || abs(dtaudx) < 10.0*ck_Hp || x_current > -8.3){
      tc = false;
    }
    else{idx_tc++;}
  }

  return idx_tc;
}

//===================================
// After integrating the perturbation
// compute the source function
//===================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=====================================================================
  // Make the x and k arrays to evaluate over and use to make the splines
  //=====================================================================
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, n_k);
  Vector k_array = exp(log_k_array);

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================
      // Compute the source functions
      //=============================
      // Fetch all the things we need
      const double c  = Constants.c;

      double Hp       = cosmo->Hp_of_x(x);
      double dHpdx    = cosmo->dHpdx_of_x(x);
      double ddHpddx  = cosmo->ddHpddx_of_x(x);

      double tau      = rec->tau_of_x(x);
      double dtaudx   = rec->dtaudx_of_x(x);
      double ddtauddx = rec->ddtauddx_of_x(x);
      double g        = rec->g_tilde_of_x(x);
      double dgdx     = rec->dgdx_tilde_of_x(x);
      double ddgddx   = rec->ddgddx_tilde_of_x(x);

      double Psi         = get_Psi(x, k);
      double dPsidx      = get_dPsidx(x, k);
      double Phi         = get_Phi(x, k);
      double dPhidx      = get_dPhidx(x, k);
      double v_b         = get_v_b(x, k);
      double dv_bdx      = get_dv_bdx(x, k);

      double Theta0      = get_Theta(x, k, 0);
      double Theta1      = get_Theta(x, k, 1);
      double Theta2      = get_Theta(x, k, 2);
      double Theta3      = get_Theta(x, k ,3);
      double dTheta1dx   = get_dThetadx(x, k, 1);
      double dTheta2dx   = get_dThetadx(x, k, 2);
      double dTheta3dx   = get_dThetadx(x, k, 3);
      
      double ck_Hp    = c*k/Hp;
      // Found analytically
      double ddTheta2ddx = 2.0/5.0*ck_Hp*(dTheta1dx - dHpdx/Hp*Theta1) - 3.0/5.0*ck_Hp*(dTheta3dx - dHpdx/Hp*Theta3) + 9.0/10.0*(ddtauddx*Theta2 + dtaudx*dTheta2dx);
      // Split up for readability
      double source_function_1 = g*(Psi + Theta0 + 1.0/4.0*Theta2) + exp(-tau)*(dPsidx - dPhidx) - 1.0/(c*k)*(dHpdx*g*v_b + Hp*dgdx*v_b + Hp*g*dv_bdx);
      double source_function_2 = 3.0*Hp/(4.0*pow(c*k, 2))*(g*(pow(dHpdx, 2)*Theta2 + 3.0*dHpdx*Hp*dTheta2dx + Hp*ddHpddx*Theta2 + pow(Hp, 2)*ddTheta2ddx) + dgdx*(3.0*dHpdx*Hp*Theta2 + 2.0*pow(Hp, 2)*dTheta2dx) + pow(Hp, 2)*ddgddx*Theta2);

      double source_function = source_function_1 + source_function_2;

      // Temperature source
      ST_array[index] = source_function;
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");

  Utils::EndTiming("source");
}

//=============================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//=============================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //===========================================================
  // Compute where in the y / dydx array each component belongs
  //===========================================================

  // The different quantities in the y array
  const double &delta_cdm      =  y[Constants.ind_deltacdm_tc];
  const double &delta_b        =  y[Constants.ind_deltab_tc];
  const double &v_cdm          =  y[Constants.ind_vcdm_tc];
  const double &v_b            =  y[Constants.ind_vb_tc];
  const double &Phi            =  y[Constants.ind_Phi_tc];
  const double *Theta          = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx         =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx           =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx             =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx               =  dydx[Constants.ind_vb_tc];
  double &dPhidx               =  dydx[Constants.ind_Phi_tc];
  double *dThetadx             = &dydx[Constants.ind_start_theta_tc];

  //================================================
  // fill in the expressions for all the derivatives
  //================================================

  // Constants, cosmo and recombination values
  const double c = Constants.c;

  double H0        = cosmo->get_H0();
  double Hp        = cosmo->Hp_of_x(x);
  double dHpdx     = cosmo->dHpdx_of_x(x);
  double OmegaR    = cosmo->get_OmegaR(x);
  double OmegaCDM0 = cosmo->get_OmegaCDM(0);
  double OmegaB0   = cosmo->get_OmegaB(0);
  double OmegaR0   = cosmo->get_OmegaR(0);

  double dtaudx    = rec->dtaudx_of_x(x);
  double ddtauddx  = rec->ddtauddx_of_x(x);
  double R         = rec->get_R(x);

  // Shorthand quantities
  double ck_Hp  = c*k/Hp;
  double Theta2 = -20.0/(45.0*dtaudx)*ck_Hp*Theta[1];
  double Psi    = -Phi - 12.*pow(H0/(c*k*exp(x)), 2)*OmegaR0*Theta2;
  

  // ODE setup
  dPhidx       = Psi - 1.0/3.0*pow(ck_Hp, 2)*Phi + 1.0/2.0*pow(H0/Hp, 2)*(OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)*Theta[0]);
  ddelta_cdmdx = ck_Hp*v_cdm - 3.0*dPhidx;
  ddelta_bdx   = ck_Hp*v_b - 3.0*dPhidx;
  dv_cdmdx     = -v_cdm - ck_Hp*Psi;

  dThetadx[0]  = -ck_Hp*Theta[1] - dPhidx;
  double q     = (-((1.0 - R)*dtaudx + (1.0 + R)*ddtauddx) * (3.0*Theta[1] + v_b) - ck_Hp*Psi + (1.0 - dHpdx/Hp)*ck_Hp*(-Theta[0] + 2.0*Theta2) - ck_Hp*dThetadx[0]) / ((1.0 + R)*dtaudx + dHpdx/Hp - 1.0);
  dv_bdx       = 1.0/(1.0 + R) * (-v_b - ck_Hp*Psi + R*(q + ck_Hp*(-Theta[0] + 2.0*Theta2 - Psi)));

  dThetadx[1]  = 1.0/3.0*(q - dv_bdx);

  return GSL_SUCCESS;
}

//====================================
// The right hand side of the full ODE
//====================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  //===========================================================
  // Compute where in the y / dydx array each component belongs
  //===========================================================
  
  // Index and number of the different quantities
  const int n_ell_theta   = Constants.n_ell_theta;

  // The different quantities in the y array
  const double &delta_cdm =  y[Constants.ind_deltacdm];
  const double &delta_b   =  y[Constants.ind_deltab];
  const double &v_cdm     =  y[Constants.ind_vcdm];
  const double &v_b       =  y[Constants.ind_vb];
  const double &Phi       =  y[Constants.ind_Phi];
  const double *Theta     = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx   =  dydx[Constants.ind_deltab];
  double &dv_cdmdx     =  dydx[Constants.ind_vcdm];
  double &dv_bdx       =  dydx[Constants.ind_vb];
  double &dPhidx       =  dydx[Constants.ind_Phi];
  double *dThetadx     = &dydx[Constants.ind_start_theta];

  // Constants, cosmo and recombination values
  const double c = Constants.c;

  double H0        = cosmo->get_H0();
  double Hp        = cosmo->Hp_of_x(x);
  double dHpdx     = cosmo->dHpdx_of_x(x);
  double eta       = cosmo->eta_of_x(x);
  double OmegaCDM0 = cosmo->get_OmegaCDM();
  double OmegaB0   = cosmo->get_OmegaB();
  double OmegaR0   = cosmo->get_OmegaR();

  double dtaudx    = rec->dtaudx_of_x(x);
  double ddtauddx  = rec->ddtauddx_of_x(x);
  double R         = rec->get_R(x);

  // Shorthand quantities
  double ck_Hp = c*k/Hp;
  double Psi   = -Phi - 12.0*pow(H0/(c*k*exp(x)), 2)*OmegaR0*Theta[2];

  //================================================
  // Fill in the expressions for all the derivatives
  //================================================

  // Scalar quantities (Phi, delta, v, ...)
  dPhidx       = Psi - 1.0/3.0*pow(ck_Hp, 2)*Phi + 1.0/2.0*pow(H0/Hp, 2)*(OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)*Theta[0]);
  ddelta_cdmdx = ck_Hp*v_cdm - 3.0*dPhidx;
  ddelta_bdx   = ck_Hp*v_b - 3.0*dPhidx;
  dv_cdmdx     = -v_cdm - ck_Hp*Psi;
  dv_bdx       = -v_b - ck_Hp*Psi + dtaudx*R*(3.0*Theta[1] + v_b);

  // Photon multipoles (Theta_ell)
  dThetadx[0] = -ck_Hp*Theta[1] - dPhidx;
  dThetadx[1] = 1.0/3.0*ck_Hp*Theta[0] - 2.0/3.0*ck_Hp*Theta[2] + 1.0/3.0*ck_Hp*Psi + dtaudx*(Theta[1] + 1.0/3.0*v_b);
  dThetadx[2] = 2.0/5.0*ck_Hp*Theta[1] - 3.0/5.0*ck_Hp*Theta[3] + 9.0/10.0*dtaudx*Theta[2];
  // For 2 < ell < l_max
  for (int ell = 3; ell < n_ell_theta - 1; ell++){
    dThetadx[ell] = ell*ck_Hp/(2.0*ell + 1.0)*Theta[ell - 1] - ck_Hp*(ell + 1.0)/(2.0*ell + 1.0)*Theta[ell + 1] + dtaudx*Theta[ell];
  }
  // For ell = l_max
  int l = n_ell_theta - 1;
  dThetadx[l] = ck_Hp*Theta[l - 1] - c*(l + 1.0)/(Hp*eta)*Theta[l] + dtaudx*Theta[l];

  return GSL_SUCCESS;
}

//============
// Get methods
//============

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}

double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}

double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}

double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}

double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}

double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}

double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}

double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}

double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}

double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}

double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}

double Perturbations::get_dThetadx(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_x(x,k);
}

double Perturbations::get_eta_k(const double x, const double k) const{
  return cosmo->eta_of_x(x)*k;
}

//=======================================
// Print some useful info about the class
//=======================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:        "  << x_start                      << "\n";
  std::cout << "x_end:           " << x_end                        << "\n";
  std::cout << "n_x:             " << n_x                          << "\n";
  std::cout << "k_min (1/Mpc):   " << k_min * Constants.Mpc        << "\n";
  std::cout << "k_max (1/Mpc):   " << k_max * Constants.Mpc        << "\n";
  std::cout << "n_k:             " << n_k                          << "\n";
  std::cout << "                 " << " "                          << "\n";
  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:    " << Constants.ind_deltacdm       << "\n";
  std::cout << "ind_deltab:      " << Constants.ind_deltab         << "\n";
  std::cout << "ind_v_cdm:       " << Constants.ind_vcdm           << "\n";
  std::cout << "ind_v_b:         " << Constants.ind_vb             << "\n";
  std::cout << "ind_Phi:         " << Constants.ind_Phi            << "\n";
  std::cout << "ind_start_theta: " << Constants.ind_start_theta    << "\n";
  std::cout << "n_ell_theta:     " << Constants.n_ell_theta        << "\n";
  std::cout << "n_ell_tot_full:  " << Constants.n_ell_tot_full     << "\n";
  std::cout << "                 " << " "                          << "\n";
  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:    " << Constants.ind_deltacdm_tc    << "\n";
  std::cout << "ind_deltab:      " << Constants.ind_deltab_tc      << "\n";
  std::cout << "ind_v_cdm:       " << Constants.ind_vcdm_tc        << "\n";
  std::cout << "ind_v_b:         " << Constants.ind_vb_tc          << "\n";
  std::cout << "ind_Phi:         " << Constants.ind_Phi_tc         << "\n";
  std::cout << "ind_start_theta: " << Constants.ind_start_theta_tc << "\n";
  std::cout << "n_ell_theta:     " << Constants.n_ell_theta_tc     << "\n";
  std::cout << "n_ell_tot_tc:    " << Constants.n_ell_tot_tc       << "\n";
  std::cout << std::endl;
}

//===================================================
// Output some results to file for a given value of k
//===================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, 0.0, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                                              << " ";
    fp << get_Theta(x, k, 0)                             << " ";
    fp << get_Theta(x, k, 1)                             << " ";
    fp << get_Theta(x, k, 2)                             << " ";
    fp << get_delta_cdm(x, k)                            << " ";
    fp << get_delta_b(x, k)                              << " ";
    fp << get_v_cdm(x, k)                                << " ";
    fp << get_v_b(x, k)                                  << " ";
    fp << get_Phi(x, k)                                  << " ";
    fp << get_Psi(x, k)                                  << " ";
    fp << get_Source_T(x, k)                             << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(5,   arg)    << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(50,  arg)    << " ";
    fp << get_Source_T(x, k) * Utils::j_ell(500, arg)    << " ";
    fp << get_eta_k(x, k)                                << " ";
    fp << rec->Xe_of_x(x)                                << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}