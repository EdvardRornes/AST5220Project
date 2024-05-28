#ifndef _POWERSPECTRUM_HEADER
#define _POWERSPECTRUM_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <functional>
#include <utility> 
#include <fstream> 
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class PowerSpectrum {
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
    Perturbations *pert        = nullptr;

    Vector make_linspace_from_dx(double min, double max, double delta, bool logarithm = false);
    double trapezoidal_integral(double dx, const std::vector<double>& y_array);

    // Parameters defining the primordial power-spectrum
    double A_s        = 2.1e-9;
    double n_s        = 0.965;
    double kpivot_mpc = 0.05;

    // The k-values we compute Theta_ell(k) etc. for
    const int n_k      = 100;
    const double k_min = Constants.k_min;
    const double k_max = Constants.k_max;

    // Sampling values
    const double n_bessel       = 25.0;  // Samplings per oscillation for the Bessel function
    const double n_k_theta_LOS  = 32.0;  // Samplings per oscillation of theta_ell
    const double n_x_LOS        = 100.0; // Samplings per oscillation of the integrand in the line of sight (LOS) integral
    const double n_k_PS         = 32.0;  // Sampling per oscillation when integrating across k in the power spectrum (PS) integral

    const double eta0   = cosmo->eta_of_x(0);
    
    const double x_start     = Constants.x_start;
    const double x_end       = Constants.x_end;
    const double x_start_LOS = -8.0;     // Time before recombination from which we integrate the line of sight integral
    const double x_end_LOS   = 0.0;
    
    // The ells's we will compute Theta_ell and Cell for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};
    // Test l = {2, 5, 10, 20, 40, 100, 200, 500, 1000}
    Vector test_ell_index{0, 3, 7, 10, 13, 19, 24, 32, 42};
   
    //=====================================================================
    // [1] Create bessel function splines needed for the LOS integration
    //=====================================================================

    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_splines;
    
    // Generate splines of bessel-functions for each ell needed
    // to do the LOS integration
    void generate_bessel_function_splines();
    
    //=====================================================================
    // [2] Do the line of sight integration and spline the result
    //=====================================================================
    
    // Do LOS integration for all ells and all k's in the given k_array
    // and for all the source functions (temperature, polarization, ...)
    void line_of_sight_integration(Vector & k_array, std::function<double(double, double)> &source_function);
  
    // Do the line of sight integration for a single quantity
    // for all ells by providing a source_function(x,k) (can be temp, pol, ...)
    Vector2D line_of_sight_integration_single(
        Vector & k_array, 
        std::function<double(double,double)> &source_function);
    
    // Splines of the reusult of the LOS integration
    std::vector<Spline> thetaT_ell_of_k_spline;
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================
    
    // General method to solve for Cells (allowing for cross-correlations)
    // For auto spectrum (C_TT) then call with f_ell = g_ell = theta_ell
    Vector solve_for_cell(
        Vector & logk_array,
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell);

    // Splines with the power-spectra
    Spline cell_TT_spline{"cell_TT_spline"};

  public:

    // Constructors
    PowerSpectrum() = delete;
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert,
        double A_s,
        double n_s,
        double kpivot_mpc);
    
    // Do all the solving: bessel functions, LOS integration and then compute Cells
    void solve();

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k)
    double primordial_power_spectrum(const double k) const;

    // Get P(k,x) for a given x in units of (Mpc)^3
    double get_matter_power_spectrum(const double x, const double k_mpc) const;

    // Get the quantities we have computed
    double get_cell_TT(const double ell) const;

    double get_bessel_func(const int il_index, const double z) const;

    double get_thetaT_ell_of_k_spline(const int il_index, const double k) const;

    // Output Cells in units of l(l+1)/2pi (muK)^2
    void output(std::string filename) const;
    void output_matter_PS(std::string filename) const;
    void output_theta(std::string filename) const;
    void output_bessel_function(std::string filename) const;
    // void output_LOS_integrand(std::string filename) const;
};

#endif