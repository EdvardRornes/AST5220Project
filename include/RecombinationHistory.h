#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = Constants.x_end;
    
    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 23001;
  
    // Xe and x for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;
    double x_transition = 0;

    // Constants from the cosmology class
    double OmegaB; 
    double TCMB;
    double H0;
    double rho_c0;   // Critical density today
    void set_cosmo_constant();

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();

    // Solve for Sound Horizon
    void solve_for_sound_horizon();

    // Splines contained in this class
    Spline Xe_of_x_spline{"Xe"};
    Spline Xe_saha_of_x_spline{"Xe Saha"};
    Spline log_ne_of_x_spline{"ne"};
    Spline tau_of_x_spline{"tau"}; 
    Spline dtaudx_of_x_spline{"dtau"};
    Spline g_tilde_of_x_spline{"g"};
    Spline dg_tildedx_of_x_spline{"dgdx"};
    Spline sound_horizon_of_x_spline{"s"};

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;

    void events(const std::string filename) const;

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double get_c_s_of_x(double x) const;
    double sound_horizon_of_x(double x) const;
    double Xe_of_x(double x) const;
    double Xe_saha_of_x(double x) const;
    double ne_of_x(double x) const;
    double get_Yp() const;
    double get_t_of_x(double x) const;
    double get_R(double x) const;
};

#endif
