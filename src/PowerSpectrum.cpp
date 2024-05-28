#include "../include/PowerSpectrum.h"


//=============
// Constructors
//=============

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

// Function to create a linspace with desired spacing
Vector PowerSpectrum::make_linspace_from_dx(double min, double max, double delta, bool logarithm) {
  int npts = abs(((max-min)/delta));
  // Allow function to get the desired npts with log spacing
  if (logarithm){
    return Utils::linspace(log(min), log(max), npts);
  }
  else {
    return Utils::linspace(min, max, npts);
  }
}

// Rudimentary integration method
double PowerSpectrum::trapezoidal_integral(double dx, const std::vector<double>& y_array) {
    double integral = 0.0;
    // Calculate integral using the trapezoid method
    for (int i = 0; i < y_array.size() - 1; i++) {
        integral += 0.5*(y_array[i] + y_array[i + 1]);
    }
    return integral*dx;
}

//===================
// Do all the solving
//===================
void PowerSpectrum::solve(){

  // generate_bessel_function_splines
  generate_bessel_function_splines();

  // Make k array for line of sight (LOS) integration
  double dk_theta = 2.0*M_PI/(eta0*n_k_theta_LOS);
  Vector k_theta_array = make_linspace_from_dx(k_min, k_max, dk_theta);

  // Get source function from perturbations
  std::function<double(double, double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  
  // Implement line_of_sight_integration
  line_of_sight_integration(k_theta_array, source_function_T);

  // Make k array for power spectrum (PS) integration
  double dk_PS = 2.0*M_PI/(eta0*n_k_PS);
  Vector log_k_PS_array = make_linspace_from_dx(k_min, k_max, dk_PS, true);

  // solve_for_cell
  auto cell_TT = solve_for_cell(log_k_PS_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
}

//========================================================
// Generate splines of j_ell(z) needed for LOS integration
//========================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  const int size_ell = ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //==============================================
  // Compute splines for bessel functions j_ell(z)
  //==============================================

  double z_min = 0.0;
  double z_max = k_max*eta0;
  double dz = 2*M_PI/n_bessel;
  Vector z_array = make_linspace_from_dx(z_min, z_max, dz);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int l = 0; l < ells.size(); l++){
    const int ell = ells[l];
    // Vector to store j_ell
    Vector j_ell_array(z_array.size());

    for (int z = 0; z < z_array.size(); z++){
      j_ell_array[z] = Utils::j_ell(ell, z_array[z]);
    }
    j_ell_splines[l].create(z_array, j_ell_array);
  }

  Utils::EndTiming("besselspline");
}

//=================================
// Do the line of sight integration
// for a single source function
//=================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double dx = 2.0*M_PI/n_x_LOS;
  Vector x_array = make_linspace_from_dx(x_start_LOS, x_end_LOS, dx);

  #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < k_array.size(); ik++){
    if (10*ik/k_array.size() != (10*ik + 10)/k_array.size()){
      std::cout << (100*ik + 100)/k_array.size() << "% " << std::flush;
      if (ik == k_array.size() - 1){
        std::cout << std::endl;
      }
    }
    double k_value = k_array[ik];
    for (int il = 0; il < ells.size(); il++){
      double ell = ells[il];

      Vector integrand(x_array.size());
      for (int i = 0; i < x_array.size(); i++){
        integrand[i] = source_function(x_array[i], k_value) * j_ell_splines[il](k_value*(eta0 - cosmo->eta_of_x(x_array[i])));
      }
      result[il][ik] = trapezoidal_integral(dx, integrand);
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//=================================
// Do the line of sight integration
//=================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array, std::function<double(double, double)> &source_function){
  const int n_ells = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(n_ells);

  // Solve for Theta_ell(k)
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function);

  // Spline the result
  for (int il = 0; il < ells.size(); il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il]);
  }
}

//=============
// Compute Cell
//=============
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int n_ells      = ells.size();

  //====================================================
  // Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  //====================================================
  Vector result(n_ells);
  int n = log_k_array.size();
  double dk_log = (log_k_array[n-1] - log_k_array[0])/n;

  // Integrate for all ells
  for (int il = 0; il < n_ells; il++){
    double ell = ells[il];
      Vector integrand(log_k_array.size());
      for (int i = 0; i < log_k_array.size(); i++){
        double k_value = exp(log_k_array[i]);
        integrand[i] = primordial_power_spectrum(k_value) * abs(f_ell_spline[il](k_value)*g_ell_spline[il](k_value));
      }
      result[il] = 4.0*M_PI*trapezoidal_integral(dk_log, integrand);
    }
  return result;
}

//======================================
// Compute the primordial power-spectrum
//======================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s*pow(Constants.Mpc*k/kpivot_mpc, n_s - 1.0);
}

//=========================
// P(k) in units of (Mpc)^3
//=========================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  //==================================
  // Compute the matter power spectrum
  //==================================
  double c       = Constants.c;
  double OmegaM0 = cosmo->get_OmegaB() + cosmo->get_OmegaCDM();
  double H0      = cosmo->get_H0();
  double Phi     = pert->get_Phi(x, k_mpc);
  
  double P_prim_matter = 2.0*pow(M_PI, 2.0)*primordial_power_spectrum(k_mpc)/(pow(k_mpc, 3.0));
  double Delta_M = 2.0*pow(c*k_mpc, 2.0)*Phi / (3.0*OmegaM0*exp(-x)*pow(H0, 2.0));

  double pofk = pow(Delta_M, 2.0)*P_prim_matter;

  return pofk;
}

//============
// Get methods
//============
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}

double PowerSpectrum::get_bessel_func(const int il_index, const double z) const{
  return j_ell_splines[il_index](z);
}

double PowerSpectrum::get_thetaT_ell_of_k_spline(const int il_index, const double k) const{
  return thetaT_ell_of_k_spline[il_index](k);
}

//====================
// Output data to file
//====================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size() - 1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell*(ell + 1)) / (2.0*M_PI)*pow(1e6*cosmo->get_TCMB(0.0), 2.0);
    fp << ell                                 << " ";
    fp << cell_TT_spline(ell) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_matter_PS(std::string filename) const{
  std::ofstream fp(filename.c_str());
  double Mpc = Constants.Mpc;
  double h = cosmo -> get_h();
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000);
  Vector k_array = exp(log_k_array);

  auto print_data = [&] (const double k) {
    fp << k*Mpc/h << " ";
    fp << get_matter_power_spectrum(0.0, k)*pow(h/Mpc, 3.) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_theta(std::string filename) const{
  // Output of theta l for l={6, 100, 200, 500, 1000}
  std::ofstream fp(filename.c_str());
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);
  Vector log_k_array = Utils::linspace(log_k_min, log_k_max, 10000); // Linearly spaced logarithmic values
  Vector k_array = exp(log_k_array);
  double c = Constants.c;
  double H0 = cosmo -> get_H0();

  auto print_data = [&] (const double k){
    fp << k*eta0 << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[0], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[1], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[2], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[3], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[4], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[5], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[6], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[7], k) << " ";
    fp << get_thetaT_ell_of_k_spline(test_ell_index[8], k) << " ";
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

void PowerSpectrum::output_bessel_function(std::string filename) const{
  // Output bessel functions
  std::ofstream fp(filename.c_str());
  Vector z_array = Utils::linspace(0.0, 1e3, 1e4);
  auto print_data = [&] (const double z) {
    fp << z << " ";
    fp << get_bessel_func(test_ell_index[0], z) << " ";
    fp << get_bessel_func(test_ell_index[1], z) << " ";
    fp << get_bessel_func(test_ell_index[2], z) << " ";
    fp << get_bessel_func(test_ell_index[3], z) << " ";
    fp << get_bessel_func(test_ell_index[4], z) << " ";
    fp << get_bessel_func(test_ell_index[5], z) << " ";
    fp << get_bessel_func(test_ell_index[6], z) << " ";
    fp << get_bessel_func(test_ell_index[7], z) << " ";
    fp << get_bessel_func(test_ell_index[8], z) << " ";
    fp << "\n";
  };
  std::for_each(z_array.begin(), z_array.end(), print_data);
}
