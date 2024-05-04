#include "../include/PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

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
  double npts = (int)abs(((max-min)/delta));
  if(logarithm){
    return Utils::linspace(log(min), log(max), npts);
  }
  else{
    return Utils::linspace(min, max, npts);
  }
}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  // Make k array for line of sight (LOS) integration
  double log_k_min = log(k_min);
  double log_k_max = log(k_max);

  double dk_LOS = 2.0*M_PI/(eta0*n_k_theta_LOS);

  Vector k_theta_array = make_linspace_from_dx(k_min, k_max, dk_LOS);
  Vector log_k_theta_array = log(k_theta_array);

  // Make k array for power spectrum (PS) integration
  double dk_PS = 2.0*M_PI/(eta0*n_k_PS);

  Vector log_k_PS_array = make_linspace_from_dx(k_min, k_max, dk_PS, true);
  Vector k_PS_array = exp(log_k_PS_array);

  //========================================================================= 
  // generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_theta_array);

  //=========================================================================
  // solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_theta_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");

  const int size_ell = (int)ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  double z_min = 0.0;
  double z_max = k_max*eta0;
  double dz = 2*M_PI/n_bessel;
  Vector z_array = make_linspace_from_dx(z_min, z_max, dz);

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];

    Vector j_ell_array(z_array.size());

    for (size_t i = 0; i < z_array.size(); i++){
      j_ell_array[i] = Utils::j_ell(ell, z_array[i]);
    }
    j_ell_splines[i].create(z_array, j_ell_array);
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double dx = 2.0*M_PI/n_x_LOS;

  Vector x_array = make_linspace_from_dx(x_start_LOS, x_end_LOS, dx);

  for(size_t ik = 0; ik < k_array.size(); ik++){
    if (10.0*ik/k_array.size() != (10.0*ik + 10)/k_array.size()){
      std::cout << (100*ik + 100)/k_array.size() << "% " << std::flush;
      if (ik == k_array.size() - 1){
        std::cout << std::endl;
      }
    }
    double k_value = k_array[ik];
    for (size_t il = 0; il < ells.size(); il++){
      double ell = ells[il];

      Vector integrand(x_array.size());
      for (size_t i = 0; i < x_array.size(); i++){
        integrand[i] = source_function(x_array[i], k_value) * j_ell_splines[il](k_value*(eta0 - cosmo->eta_of_x(x_array[i])));
      }
      double int_value = 0;
      size_t N = integrand.size();
      for (size_t i = 1; i < N - 1; i++){
        int_value += integrand[i];
      }
      int_value += (integrand[0] + integrand[N-1])/2.0;
      result[il][ik] = int_value * dx;
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int n_ells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(n_ells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================
  // Vector2D thetaT_ell_of_k_spline = line_of_sight_integration_single(k_array, pert->get_Source_T(x,k))

  // Make a function returning the source function
  // std::function<double(double,double)> source_function_T = [&](double x, double k){
  //   return pert->get_Source_T(x,k);
  // };

  // Do the line of sight integration
  // Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // ...
    // ...
    // ...
    // ...

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  Vector result;

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

