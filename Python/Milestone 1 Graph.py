import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import UnivariateSpline
from scipy.stats import norm
from scipy.optimize import curve_fit
from tabulate import tabulate

# Constants
c              = 3*10**8 
MPC_to_meter   = 1/3.086e22
Hkm_per_parsec = 10**5*MPC_to_meter
Gigayear       = 10**9*365*24*3600 # Gigayear to seconds
H0             = 0.67*10**5*MPC_to_meter

# I do not know why the following lines being the way they are is necessary but I could not run the script otherwise
# Get the current directory of the script
current_directory            = os.path.dirname(os.path.abspath(__file__))

# Construct the paths to the files using the current directory
supernova_data               = os.path.join(current_directory, '..', 'data', 'supernovadata.txt')
cosmology                    = os.path.join(current_directory, '..', 'data', 'cosmology.txt')
numerical_supernova          = os.path.join(current_directory, '..', 'data', 'results_supernovafitting.txt')

# Extract from data
data_cosmology               = np.genfromtxt(cosmology)
data_numerical_supernova     = np.genfromtxt(numerical_supernova, delimiter=None, skip_header=0)
data_supernova_data          = np.genfromtxt(supernova_data)

# Extract columns from cosmology
x_values                     = data_cosmology[:, 0]
eta_values                   = data_cosmology[:, 1]
t_values                     = data_cosmology[:, 2]
Hp_values                    = data_cosmology[:, 3]
dHpdx_values                 = data_cosmology[:, 4]
ddHpddx_values               = data_cosmology[:, 5]
cosmology_OmegaB_values      = data_cosmology[:, 6]
cosmology_OmegaCDM_values    = data_cosmology[:, 7]
cosmology_OmegaLambda_values = data_cosmology[:, 8]
cosmology_OmegaR_values      = data_cosmology[:, 9]
cosmology_OmegaNu_values     = data_cosmology[:, 10]
cosmology_OmegaK_values      = data_cosmology[:, 11]
luminosity_distance_values   = data_cosmology[:, 12]
comoving_distance_values     = data_cosmology[:, 13]
r_values                     = data_cosmology[:, 14]
cosmology_OmegaM_values      = cosmology_OmegaB_values + cosmology_OmegaCDM_values  # Non-relativistic particles
cosmology_OmegaRel_values    = cosmology_OmegaR_values + cosmology_OmegaNu_values   # Relativistic particles
a_of_x                       = np.exp(x_values)                                     # Use scale factor instead of x
redshift                     = np.exp(-x_values)-1

# Extract columns from supernova
chi_squared_values           = data_numerical_supernova[299:, 0]
h_values                     = data_numerical_supernova[299:, 1]
supernova_OmegaM_values      = data_numerical_supernova[299:, 2]
supernova_OmegaK_values      = data_numerical_supernova[299:, 3]
supernova_OmegaDE_values     = 1 - (supernova_OmegaM_values + supernova_OmegaK_values)  # Dark energy

# Extract observational supernova data
redshift_data                = data_supernova_data[:, 0]
luminosity_distance_data     = data_supernova_data[:, 1]
error_values                 = data_supernova_data[:, 2]

# Index finder function
def find_index(array, type, threshold, start, end):
    for i in range(start, end):
        if type == 0:
            if array[i] == threshold:
                output = x_values[i]
                return output
            else:
                output = None
        elif type == 1:
            if array[i] >= threshold:
                output = x_values[i]
                return output
            else:
                output = None
        elif type == -1:
            if np.abs(array[i]) <= threshold:
                output = x_values[i]
                return output
            else:
                output = None
        else:
            print('Not valid condition type')
            break

# Find the first value where matter and radiation are approximately equal
x_rad_is_mat = find_index(cosmology_OmegaM_values-cosmology_OmegaRel_values, -1, 5e-5, 0, len(x_values))

# Skip radiation domination and find the value where matter and dark energy are approximately equal.
# Notice that the threshhold is lower here to account for the steeper curves at this point.
x_mat_is_DE = find_index(cosmology_OmegaM_values - cosmology_OmegaLambda_values, -1, 5e-4, 5000, len(x_values))

# Find the value where acceleration starts. 
# Since Hp = dot(a) then dHpdx = ddot(a)/H, since H is strictly positive then dHpdx changes sign in the same way as ddot(a).
x_uni_acc = find_index(dHpdx_values, 1, 0, 0, len(x_values))

# Choose the conditions for domination to be over
arb_dom_over = 0.85
# Find the (arbitrary) first value where radiation is less than arb_dom_ver
x_rad_is_mat_arb = find_index(cosmology_OmegaRel_values, -1, arb_dom_over, 0, len(x_values))
# Skip radiation domination and find the (artibrary) first value where matter is less than arb_dom_ver
x_mat_is_DE_arb = find_index(cosmology_OmegaM_values, -1, arb_dom_over, 12000, len(x_values))


# Plot Omegai and vertical lines for various events
def Omega_of_a():
    plt.figure(figsize=(10, 6))
    plt.plot(a_of_x, cosmology_OmegaRel_values, lw = 2.5, label=r'$\Omega_\text{rad}$')
    plt.plot(a_of_x, cosmology_OmegaM_values, lw = 2.5, label=r'$\Omega_\text{M}$')
    plt.plot(a_of_x, cosmology_OmegaLambda_values, lw = 2.5, label=r'$\Omega_\Lambda$')
    plt.plot(a_of_x, cosmology_OmegaK_values, lw = 2.5, label=r'$\Omega_\text{K}$')
    plt.plot(a_of_x, cosmology_OmegaRel_values + cosmology_OmegaM_values + cosmology_OmegaLambda_values + cosmology_OmegaK_values, lw = 2.5, label = r'$\sum_i\Omega_i$')
    plt.axvline(x=np.exp(x_mat_is_DE), color='b', linestyle='--', lw = 1.5, label=r'$\Omega_{\text{rel}}(a)=\Omega_{\text{M}}(a)$')
    plt.axvline(x=np.exp(x_rad_is_mat), color='r', linestyle='--', lw = 1.5, label = r'$\Omega_{\text{M}}(a)=\Omega_\Lambda(a)$')
    plt.axvline(x=np.exp(x_uni_acc), color='g', linestyle='--', lw = 1.5, label = 'Universe accelerates')
    plt.xlim(np.exp(-18),np.exp(3))
    plt.xscale('log')
    plt.ylim(-0.05, 1.05)
    plt.xlabel(r'Scale factor $a$')
    plt.ylabel(r'$\Omega_i(x)$')
    plt.title(r'Energy density of the various $\Omega_i(a)$')
    plt.legend()
    plt.grid(True)
    plt.show()

# Test out data compared to analytical expressions for the derivatives of Hp
def derivs_of_Hp_vs_analytic():
    plt.figure(figsize=(10, 6))
    plt.axhline(y=-1, color='b', linestyle='--', lw = 2, label=r'$\frac{1}{\mathcal{H}_{\text{rel}}}\frac{d\mathcal{H}_{\text{rel}}}{dx}\approx-1$')
    plt.axhline(y=-1/2, color='g', linestyle='--', lw = 2, label=r'$\frac{1}{\mathcal{H}_{\text{M}}}\frac{d\mathcal{H}_{\text{M}}}{dx}\approx-1/2$')
    plt.axhline(y=1, color='r', linestyle='--', lw = 2, label=r'$\frac{1}{\mathcal{H}_\Lambda}\frac{d\mathcal{H}_\Lambda}{dx}\approx\frac{1}{\mathcal{H}_{\text{rel}}}\frac{d^2\mathcal{H}_{\text{rel}}}{dx^2}$' + '\n' + r'$\approx\frac{1}{\mathcal{H}_\Lambda}\frac{d^2\mathcal{H}_\Lambda}{dx^2}\approx1$')
    plt.axhline(y=1/4, color='y', linestyle='--', lw = 2, label=r'$\frac{1}{\mathcal{H}_{\text{M}}}\frac{d^2\mathcal{H}_{\text{M}}}{dx^2}\approx1/4$')
    plt.plot(a_of_x, dHpdx_values/Hp_values, lw = 3, label = r'Numerical $\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$')
    plt.plot(a_of_x, ddHpddx_values/Hp_values, lw = 3, label = r'Numerical $\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$')
    plt.xlabel(r'Scale factor $a$')
    plt.xscale('log')
    plt.xlim(np.exp(-18), np.exp(5))
    plt.ylabel(r'$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx},\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$')
    plt.title(r'Evolution of $\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$ and $\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$ compared to analytical approximations in various regimes')

    # Code to switch order of legends
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [4,5,0,1,2,3]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='center left')
    plt.grid(True)
    plt.show()


# Test out data compared to analytical expressions for eta*Hp/c in various epochs
def merge_eta_Hp_c():
    etaHpc_rel_dom        = np.ones_like(x_values)
    etaHpc_M_dom          = 2 + np.exp(x_rad_is_mat_arb - x_values / 2) - 2 * np.exp((x_rad_is_mat_arb - x_values) / 2)
    etaHpc_Lambda_dom     = (np.exp(x_rad_is_mat_arb) - 2 * np.exp(x_rad_is_mat_arb / 2) + 2 * np.exp(x_mat_is_DE_arb / 2) + np.exp(-x_mat_is_DE_arb) - np.exp(-x_values)) * np.exp(x_values)
    rad_M_eq_str          = str(x_rad_is_mat)
    rad_Lambda_eq_str     = str(x_mat_is_DE)
    rad_M_eq_arb_str      = str(x_rad_is_mat_arb)
    rad_Lambda_eq_arb_str = str(x_mat_is_DE_arb)
    fig, axs = plt.subplots(2, figsize=(10, 12))

    # Analytical approximation Omega_i = Omega_j for when epochs start and end
    axs[0].plot(x_values, eta_values * Hp_values / c, lw=3.5, label=r'$\eta\mathcal{H}/c$')
    axs[0].plot(x_values, etaHpc_rel_dom * (x_values < x_rad_is_mat), linestyle='--', lw=1.5, label=r'$\eta_{\text{rel}}\mathcal{H}_{\text{rel}}/c$')
    axs[0].plot(x_values, etaHpc_M_dom * (x_values >= x_rad_is_mat) * (x_values < x_mat_is_DE), linestyle='--', lw=1.5, label=r'$\eta_{\text{M}}\mathcal{H}_{\text{M}}/c$')
    axs[0].plot(x_values, etaHpc_Lambda_dom * (x_values >= x_mat_is_DE), linestyle='--', lw=1.5, label=r'$\eta_\Lambda\mathcal{H}_\Lambda/c$')
    axs[0].axvline(x_rad_is_mat, color=(0, 0, 0.5), lw=2.5, label=r'$\Omega_{\text{rel}}(x)=\Omega_{\text{M}}(x)$')
    axs[0].axvline(x_mat_is_DE, color=(0.5, 0, 0), lw=2.5, label=r'$\Omega_{\text{M}}(x)=\Omega_\Lambda(x)$')
    axs[0].set_xlabel(r'$x$')
    axs[0].set_xlim(-12, 1)
    axs[0].set_ylabel(r'$\eta(x)\mathcal{H}/c$')
    axs[0].set_ylim(0, 9)
    axs[0].set_title(r'$\eta\mathcal{H}/c$ compared to analytical approximations' + '\n' + r'where epochs are assumed to end abruptly once we have equality' + '\n' + r'between the density parameters i.e. $x_0=$' + rad_M_eq_str + r' and $x_1=$' + rad_Lambda_eq_str)
    axs[0].legend()
    axs[0].grid(True)

    # Analytical approximation for Omega_i < 0.8 for when epochs start and end
    axs[1].plot(x_values, eta_values * Hp_values / c, lw=3.5, label=r'$\eta\mathcal{H}/c$')
    axs[1].plot(x_values, etaHpc_rel_dom * (x_values < x_rad_is_mat_arb), linestyle='--', lw=1.5, label=r'$\eta_{\text{rel}}\mathcal{H}_{\text{rel}}/c$')
    axs[1].plot(x_values, etaHpc_M_dom * (x_values >= x_rad_is_mat_arb) * (x_values < x_mat_is_DE_arb), linestyle='--', lw=1.5, label=r'$\eta_{\text{M}}\mathcal{H}_{\text{M}}/c$')
    axs[1].plot(x_values, etaHpc_Lambda_dom * (x_values >= x_mat_is_DE_arb), linestyle='--', lw=1.5, label=r'$\eta_\Lambda\mathcal{H}_\Lambda/c$')
    axs[1].axvline(x_rad_is_mat_arb, color=(0, 0, 0.5), lw=2.5, label=r'$\Omega_{\text{rel}}(x)=\Omega_{\text{M}}(x)$')
    axs[1].axvline(x_mat_is_DE_arb, color=(0.5, 0, 0), lw=2.5, label=r'$\Omega_{\text{M}}(x)=\Omega_\Lambda(x)$')
    axs[1].set_xlabel(r'$x$')
    axs[1].set_xlim(-12, 1)
    axs[1].set_ylabel(r'$\eta(x)\mathcal{H}/c$')
    axs[1].set_ylim(0, 9)
    axs[1].set_title(r'$\eta\mathcal{H}/c$ compared to analytical approximations' + '\n' + r'where epochs are assumed to end once $\Omega_i<$'+ str(arb_dom_over) + '\n' + r'i.e. $x_0=$' + rad_M_eq_arb_str + r' and $x_1=$' + rad_Lambda_eq_arb_str)
    axs[1].legend(loc = 'upper left')
    axs[1].grid(True)

    # Adjust layout
    plt.tight_layout()
    # Adjust the vertical spacing between subplots
    plt.subplots_adjust(hspace=0.35)
    plt.show()


# Merge together plots for the evolution of Hp, t and eta as a function of x
def merge_Hp_t_eta():
    fig, axs = plt.subplots(3, 1, figsize=(10, 18))
    
    # Hp(x)
    axs[0].plot(x_values, Hp_values/Hkm_per_parsec, lw=2.5, label = r'$\mathcal{H}(x)$')
    axs[0].axvline(x_uni_acc, color='r', lw=2.5, linestyle = '--', label=r'Universe begins to accelerate $\ddot{a}\geq0$')
    axs[0].set_xlabel(r'$x$')
    axs[0].set_xlim(-12, 3)
    axs[0].set_ylabel(r'$\mathcal{H}$ [100 km/s/Mpc]')
    axs[0].set_yscale("log")
    axs[0].set_ylim(1e-1, 1.1e3)
    axs[0].set_title(r'Time evolution of the Hubble factor $\mathcal{H}(x)\equiv aH(x)$')
    axs[0].legend()
    axs[0].grid(True)
    
    # t(x)
    t_today = t_values[x_values == 0]/(Gigayear)
    t_today_str = '{:.2f}'.format(t_today[0])    
    axs[1].plot(x_values, t_values/Gigayear, lw=2.5, label=r'$t(x)/c$')
    axs[1].axhline(y=t_today, color='r', linestyle='--', lw=2.5, label=r'$t(0)=$' + t_today_str + ' Gyr')
    axs[1].set_xlabel(r'$x$')
    axs[1].set_xlim(-12, 3)
    axs[1].set_ylabel(r'$t$ [Gyr]')
    axs[1].set_yscale('log')
    axs[1].set_ylim(1e-8, 1e2)
    axs[1].set_title(r'Cosmic time $t(x)$')
    axs[1].legend()
    axs[1].grid(True)
    
    # eta(x)
    eta_today = eta_values[x_values == 0]/(c*Gigayear)
    eta_today_str = '{:.2f}'.format(eta_today[0])
    axs[2].plot(x_values, eta_values/(c*Gigayear), lw=2.5, label=r'$\eta(x)/c$')
    axs[2].axhline(y=eta_today, color='r', linestyle='--', lw=2.5, label=r'$\eta(0)/c=$' + eta_today_str + ' Gly')
    axs[2].set_xlabel(r'$x$')
    axs[2].set_xlim(-12, 3)
    axs[2].set_ylabel(r'$\eta/c$ [Gly]')
    axs[2].set_yscale('log')
    axs[2].set_ylim(5e-3, 1e2)
    axs[2].set_title(r'Conformal distance $\eta(x)/c$')
    axs[2].legend()
    axs[2].grid(True)
    
    plt.subplots_adjust(hspace=0.5)  # Adjust the vertical spacing between subplots
    
    plt.show()


# Plot the numerical luminosity distance versus luminosity distance from supernova data
def luminosity_distance():
    # Make a spline, s is how "lenient" the spline fit is. Low s means it follows the points nearly exactly
    spline = UnivariateSpline(redshift_data, luminosity_distance_data, w=1/error_values, s=50)

    # Plot luminosity distance with error bars from supernova.txt
    plt.figure(figsize=(10, 6))
    plt.errorbar(redshift_data, luminosity_distance_data/redshift_data, yerr=error_values/redshift_data, fmt='o', markersize=4, label='Data From Supernova')

    # Plot theoretical prediction, division by 0 removed
    plt.plot(redshift[redshift != 0], luminosity_distance_values[redshift != 0] / redshift[redshift != 0] * MPC_to_meter / 10**3, lw = 2.5, label='Theoretical Prediction')

    # Plot spline fit
    x_fit = np.linspace(min(x_values), max(x_values), 1000)
    plt.plot(x_fit, spline(x_fit)/x_fit, label='Spline Fit', lw = 1.5)

    plt.xlabel(r'Redshift $z$')
    plt.ylabel(r'$D_L/z$ [Gpc]')
    plt.xlim(8e-3, 1.7)
    plt.ylim(3, 8)
    plt.xscale('log')
    plt.title(r'Luminosity distance $D_L$ as a function of redshift $z$')
    plt.legend()
    plt.grid(True)
    plt.show()


# Creates and plots a scatterplot whilst extracting accepted samples
def read_and_scatterplot_data():
    try:
        # Locate the minima of chi
        chi_min = np.min(chi_squared_values)

        # Calculate confidence regions
        sigma1_mask = chi_squared_values - chi_min < 3.53
        sigma2_mask = chi_squared_values - chi_min < 8.02

        # Accepted samples
        accepted_samples_DE1 = supernova_OmegaDE_values[sigma1_mask]
        accepted_samples_DE2 = supernova_OmegaDE_values[sigma2_mask]
        accepted_samples_M1  = supernova_OmegaM_values[sigma1_mask]
        accepted_samples_M2  = supernova_OmegaM_values[sigma2_mask]
        accepted_samples_h1  = h_values[sigma1_mask]
        accepted_samples_h2  = h_values[sigma2_mask]

        # Concatenate accepted samples for posterior PDFs
        DE_accepted_samples  = np.concatenate((accepted_samples_DE1, accepted_samples_DE2))
        M_accepted_samples   = np.concatenate((accepted_samples_M1, accepted_samples_M2))
        h_accepted_samples   = np.concatenate((accepted_samples_h1, accepted_samples_h2))

        # Plotting scatter plot with confidence regions
        plt.figure(figsize=(10, 6))
        plt.scatter(accepted_samples_M2, accepted_samples_DE2, c='blue', s = 15, label=r'$2\sigma$ Constraint', rasterized = True)
        plt.scatter(accepted_samples_M1, accepted_samples_DE1, c='orange', s = 15, label=r'$1\sigma$ Constraint', rasterized = True)

        # Flat universe
        hori = np.linspace(0,1,1000)
        def f(d):
            return 1-d
        plt.plot(hori, f(hori), linestyle = '--', color = 'k', label = 'Flat universe')
        # Horizontal dashed line at y=0
        plt.axhline(1, color='r', linestyle='--')

        plt.xlabel(r'$\Omega_{\text{M}}$')
        plt.xlim(0.0, 0.7)
        plt.ylabel(r'$\Omega_{\text{DE}}$')
        plt.ylim(0.0, 1.2)

        # Code to switch order of legends
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1,0,2]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

        plt.title('Degeneracy between Dark Energy and Matter Density Parameters in Cosmology.')
        plt.show()
        
        # Returns the accepted samples
        return DE_accepted_samples, M_accepted_samples, h_accepted_samples
        
    except Exception as e:
        print("Error:", e)


# Functions to plot PDF histogram of h
def plot_histogram_h(data, bins=50):
    plt.figure(figsize=(10, 6))
    n, bins, patches = plt.hist(data, bins=bins, range=(0.675, 0.725), alpha = 0.7, density=True, color='g')

    # Fit a Gaussian distribution
    mu, std = norm.fit(data)
    mu_str = '{:.3f}'.format(mu)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 1000)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'r', linewidth=2, label = 'Gaussian fit')

    # Line at center of Gaussian
    plt.axvline(mu, color = 'b', linestyle = '--', label = r'$H_0=$' + mu_str)

    plt.xlabel(r'$H_0$ [100 km/s/Mpc]')
    plt.ylabel('Probability Density')
    plt.title(r'Probability Distribution Function of $H_0$')
    plt.legend()
    plt.grid(True)
    plt.show()


# Make a table for the values of the various time parameters in the various relevant cosmic events
def make_table():
    # Rel-Matter equality
    z_rad_is_mat            = redshift[x_values == x_rad_is_mat]
    t_rad_is_mat            = t_values[x_values == x_rad_is_mat]/Gigayear
    eta_rad_is_mat          = eta_values[x_values == x_rad_is_mat]/(c*Gigayear)
    Omega_M_rad_is_mat      = cosmology_OmegaM_values[x_values == x_rad_is_mat]
    Omega_Lambda_rad_is_mat = cosmology_OmegaLambda_values[x_values == x_rad_is_mat]
    Omega_Rel_rad_is_mat    = cosmology_OmegaRel_values[x_values == x_rad_is_mat]

    # Matter-DE equality
    z_mat_is_DE             = redshift[x_values == x_mat_is_DE]
    t_mat_is_DE             = t_values[x_values == x_mat_is_DE]/Gigayear
    eta_mat_is_DE           = eta_values[x_values == x_mat_is_DE]/(c*Gigayear)
    Omega_M_mat_is_DE       = cosmology_OmegaM_values[x_values == x_mat_is_DE]
    Omega_Lambda_mat_is_DE  = cosmology_OmegaLambda_values[x_values == x_mat_is_DE]
    Omega_Rel_mat_is_DE     = cosmology_OmegaRel_values[x_values == x_mat_is_DE]

    # Universe accelerates
    z_uni_acc               = redshift[x_values == x_uni_acc]
    t_uni_acc               = t_values[x_values == x_uni_acc]/Gigayear
    eta_uni_acc             = eta_values[x_values == x_uni_acc]/(c*Gigayear)
    Omega_M_uni_acc         = cosmology_OmegaM_values[x_values == x_uni_acc]
    Omega_Lambda_uni_acc    = cosmology_OmegaLambda_values[x_values == x_uni_acc]
    Omega_rel_uni_acc       = cosmology_OmegaRel_values[x_values == x_uni_acc]

    # Today
    x_today                 = 0.0
    z_today                 = redshift[x_values == x_today]
    t_today                 = t_values[x_values == x_today]/Gigayear
    eta_today               = eta_values[x_values == x_today]/(c*Gigayear)
    Omega_M_today           = cosmology_OmegaM_values[x_values == x_today]
    Omega_Lambda_today      = cosmology_OmegaLambda_values[x_values == x_today]
    Omega_rel_today         = cosmology_OmegaRel_values[x_values == x_today]

    table_data = [
        ['x'           , x_rad_is_mat           , x_mat_is_DE           , x_uni_acc             , x_today           ],
        ['z(x)'        , z_rad_is_mat           , z_mat_is_DE           , z_uni_acc             , z_today           ],
        ['t(x)'        , t_rad_is_mat           , t_mat_is_DE           , t_uni_acc             , t_today           ],
        ['eta(x)/c'    , eta_rad_is_mat         , eta_mat_is_DE         , eta_uni_acc           , eta_today         ],
        ['Omega_M'     , Omega_M_rad_is_mat     , Omega_M_mat_is_DE     , Omega_M_uni_acc       , Omega_M_today     ],
        ['Omega_Lambda', Omega_Lambda_rad_is_mat, Omega_Lambda_mat_is_DE, Omega_Lambda_uni_acc  , Omega_Lambda_today],
        ['Omega_Rel'   , Omega_Rel_rad_is_mat   , Omega_Rel_mat_is_DE   , Omega_rel_uni_acc     , Omega_rel_today   ]
    ]


    # Create the table
    table = tabulate(table_data, headers=["Event", "Rel-Matter Equality", "Matter-DE Equality", "Universe Accelerates", "Today"], tablefmt="simple_grid")

    print(table)
    

# Main function
def main():
    # Pick out the accepted samples whilst also plotting the scatterplot
    # Note that we only need h_accepted_samples for what is asked of in the project, but perhaps for future reference I added the others
    # DE_accepted_samples, M_accepted_samples, h_accepted_samples = read_and_scatterplot_data()
    
    # # Plot histograms
    # plot_histogram_h(h_accepted_samples)

    # # Plot various parameter combinations from cosmology file
    # luminosity_distance()
    # derivs_of_Hp_vs_analytic()
    # Omega_of_a()
    # make_table()
    # merge_Hp_t_eta()
    merge_eta_Hp_c()

main()
