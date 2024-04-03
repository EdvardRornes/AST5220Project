#include "../include/RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp):
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

// Physical constants
const double k_b         = Constants.k_b;
const double hbar        = Constants.hbar;
const double c           = Constants.c;
const double G           = Constants.G;
const double epsilon_0   = Constants.epsilon_0;
const double m_e         = Constants.m_e;
const double m_H         = Constants.m_H;
const double H0_over_h   = Constants.H0_over_h;
const double sigma_T     = Constants.sigma_T;
const double lambda_2s1s = Constants.lambda_2s1s;
const double alpha       = m_e * sqrt(3.0 / (8.0 * M_PI) * sigma_T);

void RecombinationHistory::set_cosmo_constant(){
  OmegaB = cosmo->get_OmegaB(0);
  TCMB = cosmo->get_TCMB(0);
  H0 = cosmo->get_H0();
}

const double tol = 1e-7;

void RecombinationHistory::solve(){

  // Set cosmo constants
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  solve_for_sound_horizon();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_saha(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  int last_saha_idx = 0;
  double last_saha_x = 0;
  for(int i = 0; i < npts_rec_arrays; i++){
    // Get X_e from solving the Saha equation
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit){
      saha_regime = false;
      // Store the time when we leave the Saha regime, i-1 due to us adding 1 to the index for when the check will come into play
      last_saha_idx = i-1;
      last_saha_x   = x_array[i-1];
      for (int j = last_saha_idx; j < npts_rec_arrays; j++){
        // Continue with the Saha equation for reference
        auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[j]);

        // Electron fraction and number density
        const double Xe_current = Xe_ne_data.first;
        const double ne_current = Xe_ne_data.second;
        // Remove all nonsense data where the Saha equation is not valid
        double actXe_current = Xe_current < tol ? tol: Xe_current;
        double actXe_current1 = std::isnan(actXe_current) ? tol: actXe_current;
        Xe_arr_saha[j] = actXe_current1;
      }
      break;
    }

    if(saha_regime){
      // Store the result we got from the Saha equation
      Xe_arr[i] = Xe_current;
      Xe_arr_saha[i] = Xe_current;
      ne_arr[i] = ne_current;
    }     
  }
  

  // The Peebles ODE equation
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };
    
  //=============================================================================
  // Set up IC, solve the ODE and fetch the result 
  //=============================================================================
  // New x_array for the points not solved by the Saha equation
  const int npts_peebles_array = npts_rec_arrays - last_saha_idx;
  Vector x_peebles(npts_peebles_array);

  // Fill Peebles array with Saha data
  for (int i = 0; i < npts_peebles_array; i++){
    x_peebles[i] = x_array[i + last_saha_idx];
  }

  double Xe_ini_val = Xe_arr[last_saha_idx];
  Vector Xe_ini_vec{Xe_ini_val};
  
  peebles_Xe_ode.solve(dXedx, x_peebles, Xe_ini_vec);
  auto Xe_ode = peebles_Xe_ode.get_data_by_component(0);

  // Update original Xe and ne arrays
  for(int i=last_saha_idx; i<npts_rec_arrays; i++){
    // Calculate values
    const double OmegaB0     = cosmo->get_OmegaB(0);
    const double TCMB0       = cosmo->get_TCMB(0);
    const double H0          = cosmo->get_H0();
    double Xe_temp = Xe_ode[i-last_saha_idx];
    double nb_temp = (3.0 * pow(H0, 2.0) * OmegaB0) / (8 * M_PI * G * m_H * pow(exp(x_peebles[i - last_saha_idx]), 3.0));
    double ne_temp = Xe_temp*nb_temp;

    // Fill arrays
    Xe_arr[i] = Xe_temp;
    ne_arr[i] = ne_temp;
  }

  //=============================================================================
  // Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  Vector log_ne_arr(npts_rec_arrays);

  for (int i = 0; i < npts_rec_arrays; i++){
    log_ne_arr[i] = log(ne_arr[i]);
  }

  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe");
  Xe_saha_of_x_spline.create(x_array, Xe_arr_saha, "Xe Saha");
  // Xe_saha_of_x_spline.create(x_array, saha_arr, "Xe");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);

  // Fetch cosmological parameters
  const double OmegaB0     = cosmo->get_OmegaB(0);
  const double TCMB0       = cosmo->get_TCMB(0);
  const double H0          = cosmo->get_H0();

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  double nb        = OmegaB0 * 3.0 * pow(H0, 2.0) / (8.0 * M_PI * Constants.G * m_H * pow(a, 3));
  double Tb        = TCMB0 / a;
  double prefactor = 1.0/nb * pow((k_b * m_e * Tb)/(2.0*M_PI * pow(hbar , 2)), 3.0/2.0)* exp(-epsilon_0 / (k_b * Tb));
  // Threshold to avoid overflow
  if(4.0 / prefactor < tol)
    Xe = 1.0;
  else{
    Xe = prefactor / 2.0 * (-1 + sqrt(1 + 4.0 / prefactor));
  }
  ne = nb * Xe;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);


  // Cosmological parameters
  const double OmegaB0     = cosmo->get_OmegaB(0);
  const double TCMB0       = cosmo->get_TCMB(0);
  const double H0          = cosmo->get_H0();
  double H                 = cosmo->H_of_x(x);

  double Tb                = TCMB0 / a;

  //=============================================================================
  // Write the expression for dXedx
  //=============================================================================
  double phi_2              = 0.448 * log(epsilon_0 / (Tb*k_b));
  double alpha_2            = 8.0 / sqrt(3.0 * M_PI) * c * sigma_T * sqrt(epsilon_0 / (Tb*k_b)) * phi_2;
  double beta               = alpha_2 * pow(m_e * k_b * Tb / (2.0 * M_PI * pow(hbar, 2)), 3.0/2.0) * exp(-epsilon_0 / (k_b * Tb));
  // Avoiding overflow
  double beta_2;
  if (epsilon_0 / (k_b * Tb) > 200){
    beta_2 = 0.0;
  }
  else {
    beta_2 = beta * exp(3.0 * epsilon_0 / (4.0 * k_b * Tb));
  }
  
  double n_b                = (1.0 - Yp) * (3.0 * pow(H0, 2.0) * OmegaB0) / (8 * M_PI * G * m_H * pow(a, 3.0));
  double n_H                = (1.0 - Yp) * n_b;
  double n_1s               = (1.0 - X_e) * n_H;
  double Lambda_alpha       = H * pow(3.0 * epsilon_0, 3.0) / (pow(8.0 * M_PI, 2.0) * pow(c * hbar, 3) * n_1s);
  const double Lambda_2s_1s = 8.227;
  double C_r                = (Lambda_2s_1s + Lambda_alpha) / (Lambda_2s_1s + Lambda_alpha + beta_2);
  
  double rhs = C_r / H * (beta * (1.0 - X_e) - n_H * alpha_2 * pow(X_e, 2.0));

  dXedx[0] = rhs;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("Optical Depth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 100000;
  // Make reversed array to integrate over
  Vector x_array_tau_rev = Utils::linspace(-x_end, -x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // Write the expression for dtaudx
    //=============================================================================
    // Since we are integrating backwards we must provide negative x values
    dtaudx[0] = -c * sigma_T * ne_of_x(-x) / cosmo->H_of_x(-x); 

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  Vector tau_ini_vec{0.0};

  tau_ode.solve(dtaudx, x_array_tau_rev, tau_ini_vec);

  auto tau_vec_inv = tau_ode.get_data_by_component(0);

  // Undo the reversal
  Vector x_array_tau(npts);
  Vector tau_vec(npts);
  Vector dtau_vec(npts);  // Spline does not work very well on 2nd derivatives, hence we make this one and do the first derivative.
  
  for (int i = 0; i < npts; i++){
    x_array_tau[i] = -x_array_tau_rev[npts - 1 - i];
    tau_vec[i] = -tau_vec_inv[npts - 1 - i];
  }

  for(int i = 0; i < npts; i++){
    dtau_vec[i] = -c * ne_of_x(x_array_tau[i]) * sigma_T / (cosmo->H_of_x(x_array_tau[i]));
  }

  tau_of_x_spline.create(x_array_tau, tau_vec, "tau");
  dtaudx_of_x_spline.create(x_array_tau, dtau_vec, "dtau");

  Utils::EndTiming("Optical Depth");

  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================

  Utils::StartTiming("Visibility Function");

  Vector g_tilde(npts);
  Vector dg_tildedx(npts);

  for (int i = 0; i < npts; i++){
    double const x_loc = x_array_tau[i];
    double x = x_array_tau[i];
    g_tilde[i] = -dtaudx_of_x(x_loc) * exp(-tau_of_x(x_loc));
    dg_tildedx[i] = exp(-tau_of_x(x)) * (pow(dtaudx_of_x(x), 2) - ddtauddx_of_x(x));
  }

  for(int i = 0; i < npts; i++){
    double x = x_array_tau[i];
    dg_tildedx[i] = exp(-tau_of_x(x)) * (dtaudx_of_x(x) * dtaudx_of_x(x) - ddtauddx_of_x(x));
  }

  g_tilde_of_x_spline.create(x_array_tau, g_tilde, "g");
  dg_tildedx_of_x_spline.create(x_array_tau, dg_tildedx, "dgdx");
  
  Utils::EndTiming("Visibility Function");
}

void RecombinationHistory::solve_for_sound_horizon(){

  Utils::StartTiming("Sound Horizon");
  const int npts = 100000;
  Vector x_array_sound = Utils::linspace(x_start, x_end, npts);

  ODESolver sound_ode;
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    dsdx[0] = get_c_s_of_x(x) / cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };

  double ini_val = get_c_s_of_x(x_array_sound[0]) / cosmo->Hp_of_x(x_array_sound[0]);
  Vector sound_ini_vec{ini_val};

  sound_ode.solve(dsdx, x_array_sound, sound_ini_vec);

  auto sound_vec = sound_ode.get_data_by_component(0);

  sound_horizon_of_x_spline.create(x_array_sound, sound_vec, "s");

  Utils::EndTiming("Sound Horizon");

}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::Xe_saha_of_x(double x) const{
  return Xe_saha_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::get_R(double x) const{
  return 4.0 * cosmo->get_OmegaR(0) / (3.0 * cosmo->get_OmegaB(0) * exp(x));
}

double RecombinationHistory::get_c_s_of_x(double x) const{
  return c * sqrt(get_R(x) / (3.0 * (1.0 + get_R(x))));
}

double RecombinationHistory::sound_horizon_of_x(double x) const{
  return sound_horizon_of_x_spline(x);
}


double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 1e6;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                     << " "; // 0
    fp << Xe_of_x(x)            << " "; // 1
    fp << Xe_saha_of_x(x)       << " "; // 2
    fp << ne_of_x(x)            << " "; // 3
    fp << tau_of_x(x)           << " "; // 4
    fp << dtaudx_of_x(x)        << " "; // 5
    fp << ddtauddx_of_x(x)      << " "; // 6
    fp << g_tilde_of_x(x)       << " "; // 7
    fp << dgdx_tilde_of_x(x)    << " "; // 8
    fp << ddgddx_tilde_of_x(x)  << " "; // 9
    fp << sound_horizon_of_x(x) << " "; // 10
    fp << cosmo->t_of_x(x)      << " "; // 11
    fp << cosmo->z_of_x(x)      << " "; // 12
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

// Find when recombination and last scattering occurred
// NB: I am not sure why but the "recombination_events.txt" file must be deleted before this function will create a new one. It won't overwrite the previous one.
void RecombinationHistory::events(const std::string filename) const{
  const int npts = 1e6;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // Find when Xe <= 0.1
  double x_recombination = x_array[0];
  double Xe_min = Xe_of_x(x_recombination);

  for (int i = 0; i < npts; i++){
    double Xe_temp = Xe_of_x(x_array[i]);
    if (abs(Xe_temp-0.1) < abs(Xe_min-0.1)){
      x_recombination = x_array[i];
      Xe_min = Xe_of_x(x_recombination);
    }
  }

  // Make some random guess which will of course be wrong
  double x_last_scatt = x_array[0];
  double g_max        = g_tilde_of_x(x_last_scatt);

  // Iterate through and switch out for the larger value every time we find a larger avlue
  for (int i = 0; i < npts; i++){
    double g_temp = g_tilde_of_x(x_array[i]);
    if (g_temp > g_max){
      x_last_scatt = x_array[i];
      g_max        = g_tilde_of_x(x_last_scatt);
    }
  }

  // Find tau = 1
  double x_last_scatt_tau;
  double tau_unity;
  for (int i = 0; i < npts; i++){
    double tau_temp = tau_of_x(x_array[i]);
    if (abs(tau_temp) < (1)){
      x_last_scatt_tau = x_array[i];
      tau_unity = tau_of_x(x_last_scatt_tau);
      break;
    }
  }

  // Get recombination values
  double z_recombination = cosmo->z_of_x(x_recombination);
  double t_recombination = cosmo->t_of_x(x_recombination);
  double s_recombination = sound_horizon_of_x(x_recombination);

  // Get last scattering values
  double z_last_scatt = cosmo->z_of_x(x_last_scatt);
  double t_last_scatt = cosmo->t_of_x(x_last_scatt);
  double s_last_scatt = sound_horizon_of_x(x_last_scatt);

  // Get last scattering values
  double z_last_scatt_tau = cosmo->z_of_x(x_last_scatt_tau);
  double t_last_scatt_tau = cosmo->t_of_x(x_last_scatt_tau);
  double s_last_scatt_tau = sound_horizon_of_x(x_last_scatt_tau);


  const double yr = 3600*24*365;
  const double Mpc = 3.086e22;

  // write to file
  std::ofstream fp(filename.c_str());
  fp << "Event                  " << " " << "    x       "      << " " << "    z       "      << " " << "    t [yr]  "         <<  " " << "    s [Mpc] "          << " " << "\n";
  fp << "Recombination          " << " " << x_recombination     << " " << z_recombination     << " " << t_recombination/yr     <<  " " << s_recombination/Mpc     << " " << "\n";
  fp << "Last scattering        " << " " << x_last_scatt        << " " << z_last_scatt        << " " << t_last_scatt/yr        <<  " " << s_last_scatt/Mpc        << " " << "\n";
  fp << "Last scattering w/ tau " << " " << x_last_scatt_tau    << " " << z_last_scatt_tau    << " " << t_last_scatt_tau/yr    <<  " " << s_last_scatt_tau/Mpc    << " " << "\n";

}
