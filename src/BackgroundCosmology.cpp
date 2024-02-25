#include "../include/BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB0, 
    double OmegaCDM0, 
    double OmegaK0,
    double Neff, 
    double TCMB0) :
  h(h),
  OmegaB0(OmegaB0),
  OmegaCDM0(OmegaCDM0),
  OmegaK0(OmegaK0),
  Neff(Neff), 
  TCMB0(TCMB0)
{
  //=============================================================================
  // Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================
  double pi = M_PI;
  H0 = h * Constants.H0_over_h;
  OmegaR0 = 2.0 * pow(pi, 2.0) / 30.0 * pow(Constants.k_b * TCMB0, 4.0) / (pow(Constants.hbar, 3.0) 
            * pow(Constants.c, 5.0)) * 8.0 * pi * Constants.G/(3.0 * pow(H0, 2.0));
  OmegaNu0 = Neff * 7.0 / 8.0 * pow(4.0 / 11.0, 4.0 / 3.0) * OmegaR0;
  OmegaLambda0 = 1.0-(OmegaK0 + OmegaB0 + OmegaCDM0 + OmegaR0 + OmegaNu0);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
  // Set the number of points
  const int npts = 1e4;
  // Set the range of the x-axis by using x_start and x_end from Utils.h
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    // Sets the rhs of the detadx ODE
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Set the initial condition
  double eta_ini = Constants.c / Hp_of_x(x_start);
  // Define vector eta to be used by the ODE solver
  Vector eta{eta_ini};  
  // Solve the ODE
  ODESolver eta_ode;
  eta_ode.solve(detadx, x_array, eta);
  // Pick out the relevant column
  Vector eta_array = eta_ode.get_data_by_component(0);
  // Spline the result eta_of_x_spline 
  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");    

  Utils::EndTiming("Eta");

  // The ODE for dtdx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    // Sets the rhs of the dtdx ODE
    dtdx[0] = 1.0 / H_of_x(x);
    return GSL_SUCCESS;
  };

  // Set the initial condition
  double t_ini = 1 / (2 * H_of_x(x_start));
  // Define vector t to be used by the ODE solver
  Vector t{t_ini};
  // Solve the ODE
  ODESolver t_ode;
  t_ode.solve(dtdx, x_array, t);
  // Pick out the relevant column
  Vector t_array = t_ode.get_data_by_component(0);
  // Spline the result of t_of_x_spline
  t_of_x_spline.create(x_array, t_array, "t_of_x");

}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  // Shorthand notation
  double OmegaNonRel0 = OmegaB0 + OmegaCDM0;
  double OmegaRel0    = OmegaR0 + OmegaNu0;
  // Compute H
  double H = H0 * sqrt(OmegaNonRel0 * exp(-3.0 * x) + OmegaRel0*exp(-4.0 * x)
           + OmegaK0 * exp(-2.0 * x) + OmegaLambda0);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  // Compute Hp = aH
  double Hp_of_x = exp(x) * H_of_x(x);

  return Hp_of_x;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  // Shorthand notation
  double OmegaNonRel0 = OmegaB0 + OmegaCDM0;
  double OmegaRel0    = OmegaR0 + OmegaNu0;
  // Compute the derivative of \mathcal{H} w.r.t. x. note that this formula was found analytically.
  double dHpdx_of_x = H0 * exp(-3.0*x) * (2.0 * OmegaLambda0 * exp(4.0 * x) 
                    - OmegaNonRel0 * exp(x) - 2.0 * OmegaRel0) / (2.0 * H_of_x(x) / H0);

  return dHpdx_of_x;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  // Shorthand notation
  double OmegaNonRel0 = OmegaB0 + OmegaCDM0;
  double OmegaRel0 = OmegaR0 + OmegaNu0;
  // Compute the second derivative of \mathcal{H} w.r.t. x, again found analytically.
  double ddHpddx_of_x = H0 * exp(-7.0 * x) / (4.0 * pow(H_of_x(x) / H0, 3.0)) 
                      * (4.0 * pow(OmegaLambda0, 2) * exp(8.0 * x) + 8.0 * OmegaK0 
                      * OmegaLambda0 * exp(6.0 * x) + 14.0 * OmegaNonRel0 * OmegaLambda0 
                      * exp(5.0 * x) + 24.0 * OmegaRel0 * OmegaLambda0 * exp(4.0 * x) + 2.0 
                      * OmegaNonRel0 * OmegaK0 * exp(3.0 * x) + (8.0 * OmegaRel0 * OmegaK0 + pow(OmegaNonRel0, 2.0)) 
                      * exp(2.0 * x) + 6.0 * OmegaNonRel0 * OmegaRel0 * exp(x) + 4.0 * pow(OmegaRel0, 2.0));
  

  return ddHpddx_of_x;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB0;
  // Compute OmegaB, notice that Hp_of_x is used and not H_of_x
  double OmegaB = OmegaB0 / (exp(x) * pow(Hp_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaB;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR0;
  // Compute OmegaR, notice that Hp_of_x is used and not H_of_x
  double OmegaR = OmegaR0 / (exp(2.0 * x) * pow(Hp_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaR;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu0;
  // Compute OmegaNu, notice that Hp_of_x is used and not H_of_x
  double OmegaNu = OmegaNu0 / (exp(2.0 * x) * pow(Hp_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaNu;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM0;
  // Compute OmegaCDM, notice that Hp_of_x is used and not H_of_x
  double OmegaCDM = OmegaCDM0 / (exp(x) * pow(Hp_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaCDM;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda0;
  // Compute OmegaLambda, notice that it is NOT Hp_of_x that is used, but H_of_x
  double OmegaLambda = OmegaLambda0 / (pow(H_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaLambda;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK0;
  // Compute OmegaK, notice that Hp_of_x is used and not H_of_x
  double OmegaK = OmegaK0 / (pow(Hp_of_x(x), 2.0) / pow(H0, 2.0));

  return OmegaK;
}
// Compute r for the various scenarios depending on the geometry of the universe
double BackgroundCosmology::get_r_of_x(double x) const{
  // Shorthand notation
  double shorthand = sqrt(abs(OmegaK0)) * H0 * get_comoving_distance_of_x(x) / Constants.c;
  if(x == 0.0) return r_of_x;  
  // If the universe is approximately flat
  else if(abs(OmegaK0) <= 1e-4) {
    return get_comoving_distance_of_x(x);
  }
  // If the universe is closed
  else if(OmegaK0 < -1e-4) {
    return get_comoving_distance_of_x(x) * sin(shorthand) / shorthand;
  }
  // If the universe is open
  else if(OmegaK0 > 1e-4) {
    return get_comoving_distance_of_x(x) * sinh(shorthand) / shorthand;
  }
  return r_of_x; // If the program does not detect one of these possibilities
}
// Compute the angular distance
double BackgroundCosmology::get_dA_of_x(double x) const{
  if(x == 0.0) return dA_of_x;

  double dA_of_x = exp(x) * get_r_of_x(x);

  return dA_of_x;
}

double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  // Compute the luminosity distance
  double dL_of_x = get_dA_of_x(x) * exp(-2.0 * x);

  return dL_of_x;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  // Compute the co-moving distance

  double comoving_distance_of_x = eta_of_x(0) - eta_of_x(x);

  return comoving_distance_of_x;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB0;
  return TCMB0 * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB0:         " << OmegaB0      << "\n";
  std::cout << "OmegaCDM0:       " << OmegaCDM0    << "\n";
  std::cout << "OmegaLambda0:    " << OmegaLambda0 << "\n";
  std::cout << "OmegaK0:         " << OmegaK0      << "\n";
  std::cout << "OmegaNu0:        " << OmegaNu0     << "\n";
  std::cout << "OmegaR0:         " << OmegaR0      << "\n";
  std::cout << "Neff:            " << Neff         << "\n";
  std::cout << "h:               " << h            << "\n";
  std::cout << "TCMB0:           " << TCMB0        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -18.0;
  const double x_max =  5.0;
  const int    n_pts =  23001;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                                 << " "; // 0
    fp << eta_of_x(x)                       << " "; // 1 
    fp << t_of_x(x)                         << " "; // 2
    fp << Hp_of_x(x)                        << " "; // 3
    fp << dHpdx_of_x(x)                     << " "; // 4
    fp << ddHpddx_of_x(x)                   << " "; // 5
    fp << get_OmegaB(x)                     << " "; // 6
    fp << get_OmegaCDM(x)                   << " "; // 7
    fp << get_OmegaLambda(x)                << " "; // 8
    fp << get_OmegaR(x)                     << " "; // 9
    fp << get_OmegaNu(x)                    << " "; // 10
    fp << get_OmegaK(x)                     << " "; // 11
    fp << get_luminosity_distance_of_x(x)   << " "; // 12
    fp << get_comoving_distance_of_x(x)     << " "; // 13
    fp << get_r_of_x(x)                     << " "; // 14
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

