import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.stats import norm
from tabulate import tabulate

plt.rcParams['text.usetex'] = True
plt.rcParams['axes.titlepad'] = 20 

font = {'family' : 'euclid',
        'weight' : 'bold',
        'size'   : 19}

matplotlib.rc('font', **font)

# Constants
c              = 3*10**8 
MPC_to_meter   = 1/3.086e22
Hkm_per_parsec = 10**5*MPC_to_meter
Gigayear       = 10**9*365*24*3600 # Gigayear to seconds

# Construct the paths to the files using the current directory
supernova_data               = '../data/supernovadata.txt'
cosmology                    = '../data/cosmology.txt'
numerical_supernova          = '../data/results_supernovafitting.txt'

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
x_rad_is_mat = find_index(cosmology_OmegaM_values-cosmology_OmegaRel_values, -1, 2e-5, 0, len(x_values))

# Skip radiation domination and find the value where matter and dark energy are approximately equal.
# Notice that the threshhold is lower here to account for the steeper curves at this point.
x_mat_is_DE = find_index(cosmology_OmegaM_values - cosmology_OmegaLambda_values, -1, 5e-5, 5000, len(x_values))

# Find the value where acceleration starts. 
# Since Hp = dot(a) then dHpdx = ddot(a)/H, since H is strictly positive then dHpdx changes sign in the same way as ddot(a).
x_uni_acc = find_index(dHpdx_values, 1, 0, 0, len(x_values))


# Plot Omegai and vertical lines for various events
def Omega_of_x():
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, cosmology_OmegaRel_values, lw = 2.5, label=r'$\Omega_\mathrm{rad}$')
    plt.plot(x_values, cosmology_OmegaM_values, lw = 2.5, label=r'$\Omega_\mathrm{M}$')
    plt.plot(x_values, cosmology_OmegaLambda_values, lw = 2.5, label=r'$\Omega_\Lambda$')
    plt.plot(x_values, cosmology_OmegaK_values, lw = 2.5, label=r'$\Omega_\mathrm{K}$')
    plt.plot(x_values, cosmology_OmegaRel_values + cosmology_OmegaM_values + cosmology_OmegaLambda_values + cosmology_OmegaK_values, lw = 2.5, label = r'$\sum_i\Omega_i$')
    plt.axvline(x=x_mat_is_DE, color='b', linestyle='--', lw = 2.5, label=r'$\Omega_{\mathrm{rel}}(x)=\Omega_{\mathrm{M}}(x)$')
    plt.axvline(x=x_rad_is_mat, color='r', linestyle='--', lw = 2.5, label = r'$\Omega_{\mathrm{M}}(x)=\Omega_\Lambda(x)$')
    plt.axvline(x=x_uni_acc, color='g', linestyle='--', lw = 2.5, label = 'Universe accelerates')
    plt.xlim(-18,3)
    plt.ylim(-0.05, 1.05)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Omega_i(x)$')
    plt.title(r'Energy density of the various $\Omega_i(x)$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone1/Figures/Omega_i.pdf", format='pdf')

# Test out data compared to analytical expressions for the derivatives of Hp
def derivs_of_Hp_vs_analytic():
    plt.figure(figsize=(10, 6))
    plt.axhline(y=-1, color='b', linestyle='--', lw = 2, label=r'$\frac{\mathcal{H}^\prime_{\mathrm{R}}}{\mathcal{H}_{\mathrm{R}}}=-1$')
    plt.axhline(y=-1/2, color='g', linestyle='--', lw = 2, label=r'$\frac{\mathcal{H}_{\mathrm{M}}^\prime}{\mathcal{H}_{\mathrm{M}}}=-1/2$')
    plt.axhline(y=1, color='r', linestyle='--', lw = 2, label=r'$\frac{\mathcal{H}_\Lambda^\prime}{\mathcal{H}_\Lambda}=\frac{\mathcal{H}_{\mathrm{R}}^{\prime\prime}}{\mathcal{H}_{\mathrm{R}}}=\frac{\mathcal{H}_\Lambda^{\prime\prime}}{\mathcal{H}_\Lambda}=1$')
    plt.axhline(y=1/4, color='y', linestyle='--', lw = 2, label=r'$\frac{\mathcal{H}_{\mathrm{M}}^{\prime\prime}}{\mathcal{H}_{\mathrm{M}}}=1/4$')
    plt.plot(x_values, dHpdx_values/Hp_values, lw = 3, label = r'Numerical $\frac{\mathcal{H}^{\prime}}{\mathcal{H}}$')
    plt.plot(x_values, ddHpddx_values/Hp_values, lw = 3, label = r'Numerical $\frac{\mathcal{H}^{\prime\prime}}{\mathcal{H}}$')
    plt.axvline(x_rad_is_mat, color=(0, 0, 0.5), linestyle = '-.', lw=2.5, label=r'$\Omega_{\mathrm{rel}}(x)=\Omega_{\mathrm{M}}(x)$')
    plt.axvline(x_mat_is_DE, color=(0.5, 0, 0), linestyle = '-.', lw=2.5, label=r'$\Omega_{\mathrm{M}}(x)=\Omega_\Lambda(x)$')
    plt.xlabel(r'$x$')
    plt.xlim(-18, 5)
    plt.title(r'Evolution of $\frac{\mathcal{H}^\prime}{\mathcal{H}}$ and $\frac{\mathcal{H}^{\prime\prime}}{\mathcal{H}}$')

    # Code to switch order of legends
    handles, labels = plt.gca().get_legend_handles_labels()
    order = [4,5,0,1,2,3]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='center left')
    plt.grid(True)
    plt.savefig("../Latex/Milestone1/Figures/dHpddHp_vs_anal.pdf", format='pdf')


# Test out data compared to analytical expressions for eta*Hp/c in various epochs
def eta_Hp_c():
    plt.figure(figsize=(10, 6))
    etaHpc_rel_dom        = np.ones_like(x_values)
    etaHpc_M_dom          = 2 + np.exp(x_rad_is_mat - x_values / 2) - 2 * np.exp((x_rad_is_mat - x_values) / 2)
    etaHpc_Lambda_dom     = (np.exp(x_rad_is_mat) - 2 * np.exp(x_rad_is_mat / 2) + 2 * np.exp(x_mat_is_DE / 2) + np.exp(-x_mat_is_DE) - np.exp(-x_values)) * np.exp(x_values)

    # Analytical approximation Omega_i = Omega_j for when epochs start and end
    plt.plot(x_values, eta_values * Hp_values / c, lw=3.5, label=r'$\eta\mathcal{H}/c$')
    plt.plot(x_values, etaHpc_rel_dom * (x_values < x_rad_is_mat), linestyle='--', lw=1.5, label=r'$\eta_{\mathrm{rel}}\mathcal{H}_{\mathrm{rel}}/c$')
    plt.plot(x_values, etaHpc_M_dom * (x_values >= x_rad_is_mat) * (x_values < x_mat_is_DE), linestyle='--', lw=1.5, label=r'$\eta_{\mathrm{M}}\mathcal{H}_{\mathrm{M}}/c$')
    plt.plot(x_values, etaHpc_Lambda_dom * (x_values >= x_mat_is_DE), linestyle='--', lw=1.5, label=r'$\eta_\Lambda\mathcal{H}_\Lambda/c$')
    plt.axvline(x_rad_is_mat, color=(0, 0, 0.5), lw=2.5, label=r'$\Omega_{\mathrm{rel}}(x)=\Omega_{\mathrm{M}}(x)$')
    plt.axvline(x_mat_is_DE, color=(0.5, 0, 0), lw=2.5, label=r'$\Omega_{\mathrm{M}}(x)=\Omega_\Lambda(x)$')
    plt.xlabel(r'$x$')
    plt.xlim(-12, 1)
    plt.ylabel(r'$\eta(x)\mathcal{H}/c$')
    plt.ylim(0.2, 9)
    plt.title(r'$\eta\mathcal{H}/c$ compared to analytical approximations')
    plt.legend(loc='upper left')
    plt.grid(True)

    # Adjust layout
    plt.tight_layout()
    # Adjust the vertical spacing between subplots
    plt.subplots_adjust(hspace=0.35)
    plt.savefig("../Latex/Milestone1/Figures/Eta_vs_anal_merge.pdf", format='pdf')


# Merge together plots for the evolution of Hp, t and eta as a function of x
def merge_Hp_t_eta():
    fig, axs = plt.subplots(3, 1, figsize=(10, 18))
    
    # Hp(x)
    axs[0].plot(x_values, Hp_values/Hkm_per_parsec, lw=2.5, label = r'$\mathcal{H}(x)$')
    axs[0].axvline(x_uni_acc, color='r', lw=2.5, linestyle = '--', label=r'$\ddot{a}\geq0$')
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
    axs[1].plot(x_values, t_values/Gigayear, lw=2.5, label=r'$t(x)$')
    axs[1].axhline(y=t_today, color='r', linestyle='--', lw=2.5, label=r'$t(0)=\,$' + t_today_str + ' Gyr')
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
    axs[2].axhline(y=eta_today, color='r', linestyle='--', lw=2.5, label=r'$\eta(0)/c=\,$' + eta_today_str + ' Gyr')
    axs[2].set_xlabel(r'$x$')
    axs[2].set_xlim(-12, 3)
    axs[2].set_ylabel(r'$\eta/c$ [Gyr]')
    axs[2].set_yscale('log')
    axs[2].set_ylim(5e-3, 1e2)
    axs[2].set_title(r'Conformal time $\eta(x)/c$')
    axs[2].legend()
    axs[2].grid(True)
    
    plt.subplots_adjust(hspace=0.7)  # Adjust the vertical spacing between subplots
    plt.savefig("../Latex/Milestone1/Figures/merge_Hp_t_eta_Ev.pdf", format='pdf')
    


# Plot the numerical luminosity distance versus luminosity distance from supernova data
def luminosity_distance():
    # Make a spline, s is how "lenient" the spline fit is. Low s means it follows the points nearly exactly
    spline = UnivariateSpline(redshift_data, luminosity_distance_data, w=1/error_values, s=50)

    # Plot luminosity distance with error bars from supernova.txt
    plt.figure(figsize=(10, 6))
    plt.errorbar(redshift_data, luminosity_distance_data/redshift_data, yerr=error_values/redshift_data, fmt='o', markersize=4, label='Data from supernova')

    # Plot theoretical prediction, division by 0 removed
    plt.plot(redshift[redshift != 0], luminosity_distance_values[redshift != 0] / redshift[redshift != 0] * MPC_to_meter / 10**3, lw = 2.5, label='Theoretical prediction')

    # Plot spline fit
    x_fit = np.linspace(min(x_values), max(x_values), 1000)
    plt.plot(x_fit, spline(x_fit)/x_fit, label='Spline fit to data', lw = 2.5)

    plt.xlabel(r'Redshift $z$')
    plt.ylabel(r'$D_L/z$ [Gpc]')
    plt.xlim(8e-3, 1.7)
    plt.ylim(3, 8)
    plt.xscale('log')
    plt.title(r'Luminosity distance $D_L$ as a function of redshift $z$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone1/Figures/LumiDistance.pdf", format='pdf')


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
        plt.scatter(accepted_samples_M2, accepted_samples_DE2, c='blue', s = 15, label=r'$2\sigma$ constraint', rasterized = True)
        plt.scatter(accepted_samples_M1, accepted_samples_DE1, c='orange', s = 15, label=r'$1\sigma$ constraint', rasterized = True)

        # Flat universe
        hori = np.linspace(0,1,1000)
        def f(d):
            return 1-d
        plt.plot(hori, f(hori), linestyle = '--', lw = 2.5, color = 'k', label = 'Flat universe')
        # Horizontal dashed line at y=0
        plt.axhline(1, color='r', linestyle='--', lw = 2.5)

        plt.xlabel(r'$\Omega_{\mathrm{M}}$')
        plt.xlim(0.0, 0.7)
        plt.ylabel(r'$\Omega_{\Lambda}$')
        plt.ylim(0.0, 1.2)

        # Code to switch order of legends
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1,0,2]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'lower right')

        plt.title('Degeneracy between dark energy and matter density parameters in cosmology.')
        plt.savefig("../Latex/Milestone1/Figures/ScattPlot.pdf", format='pdf')
        
        # Returns the accepted samples
        return DE_accepted_samples, M_accepted_samples, h_accepted_samples
        
    except Exception as e:
        print("Error:", e)


# Functions to plot PDF histogram of h
def plot_histogram_h(data, bins=50):
    plt.figure(figsize=(10, 6))
    n, bins, patches = plt.hist(data, bins=bins, range=(0.66, 0.74), alpha = 0.7, density=True, color='g')

    # Fit a Gaussian distribution
    mu, std = norm.fit(data)
    mu_str = '{:.3f}'.format(mu)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 1000)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'r', linewidth=2.5, label = 'Gaussian fit')

    # Line at center of Gaussian
    plt.axvline(mu, color = 'b', linestyle = '--', lw = 2.5, label = r'$H_0=\,$' + mu_str)

    # Fiducial cosmology
    plt.axvline(0.67, color = 'k', linestyle = '--', lw = 2.5, label = r'Fiducial value')

    plt.xlim(0.66, 0.74)
    plt.xlabel(r'$H_0$ [100 km/s/Mpc]')
    plt.ylabel('Probability density')
    plt.title(r'Probability distribution of $H_0$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone1/Figures/PDEh.pdf", format='pdf')


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
    DE_accepted_samples, M_accepted_samples, h_accepted_samples = read_and_scatterplot_data()
    
    # Plot histograms
    plot_histogram_h(h_accepted_samples)

    # Plot various parameter combinations from cosmology file
    luminosity_distance()
    derivs_of_Hp_vs_analytic()
    Omega_of_x()
    make_table()
    merge_Hp_t_eta()
    eta_Hp_c()
    # plt.show()  # Uncomment to show plots, otherwise it just saves them to /Latex/Milestone1/Figures

main()
