import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

# Filepaths
c_ells           = '../data/cells.txt'
planck_data_low  = '../data/planck_cell_low.txt'
planck_data_high = '../data/planck_cell_high.txt'
matter_PS        = '../data/matter_power_spectrum.txt'
transfer         = '../data/theta.txt'
cosmology        = '../data/cosmology.txt'

# Extract from data
data_c_ells = np.genfromtxt(c_ells)
data_matter_PC = np.genfromtxt(matter_PS)
data_transfer = np.genfromtxt(transfer)

# Extract columns
ell_values        = data_c_ells[:, 0]
c_ell_values      = data_c_ells[:, 1]
k_values          = data_matter_PC[:, 0]
P_k_values        = data_matter_PC[:, 1]
k_eta0_values     = data_transfer[:, 0]
theta_5_values    = data_transfer[:, 2]
theta_20_values   = data_transfer[:, 4]
theta_100_values  = data_transfer[:, 6]
theta_200_values  = data_transfer[:, 7]
theta_1000_values = data_transfer[:, 9]

# Extract from data
data_cosmology               = np.genfromtxt(cosmology)

# Extract columns from cosmology
x_values                     = data_cosmology[:, 0]
Hp_values                    = data_cosmology[:, 3]
cosmology_OmegaB_values      = data_cosmology[:, 6]
cosmology_OmegaCDM_values    = data_cosmology[:, 7]
cosmology_OmegaR_values      = data_cosmology[:, 9]
cosmology_OmegaNu_values     = data_cosmology[:, 10]
cosmology_OmegaM_values      = cosmology_OmegaB_values + cosmology_OmegaCDM_values  # Non-relativistic particles
cosmology_OmegaRel_values    = cosmology_OmegaR_values + cosmology_OmegaNu_values   # Relativistic particles


# Index finder function
def find_index(x_array, array, type, threshold, start, end):
    for i in range(start, end):
        if type == 0:
            if array[i] == threshold:
                output = x_array[i]
                return output
            else:
                output = None
        if type == 2:
            if array[i] == threshold:
                output = i
                return output
            else:
                output = None
        elif type == 1:
            if array[i] >= threshold:
                output = x_array[i]
                return output
            else:
                output = None
        elif type == -1:
            if np.abs(array[i]) <= threshold:
                output = x_array[i]
                return output
            else:
                output = None
        elif type == -2:
            if np.abs(array[i]) <= threshold:
                output = i
                return output
            else:
                output = None
        else:
            print('Not valid condition type')
            break
    
# Find the first value where matter and radiation are approximately equal
x_rad_is_mat = find_index(x_values, cosmology_OmegaM_values-cosmology_OmegaRel_values, -2, 2e-5, 0, len(x_values))

def read_data(file_path):
    data = np.loadtxt(file_path)
    x = data[:, 0]
    y = data[:, 1]
    lower_err = data[:, 2]
    upper_err = data[:, 3]
    return x, y, lower_err, upper_err

def planck_compare(filename):
    plt.figure(figsize=(10, 6))
    x, y, lower_err, upper_err = read_data(filename)
    plt.errorbar(x, y, yerr=[lower_err, upper_err], fmt='o', capsize=5)
    plt.plot(ell_values, c_ell_values, lw = 2.5, label=r'$\ell(\ell+1)C_\ell/(2\pi)$')
    plt.xlim(0, 50)
    plt.ylim(-100, 2500)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)\,[\mu\mathrm{K}]^2$')
    plt.title(r'$\ell(\ell+1)C_\ell/(2\pi)$ Compared to Planck data')
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/low_C_ell_vs_Planck.pdf", format='pdf')

def Mat_Pow():
    plt.figure(figsize=(10, 6))
    matter_PS_obs_1 = np.loadtxt('../data/matter_power_spectrum_data_1.txt')
    matter_PS_obs_2 = np.loadtxt('../data/matter_power_spectrum_data_2.txt')
    k_1 = matter_PS_obs_1[:, 0]
    P_1 = matter_PS_obs_1[:, 1]
    err_1 = matter_PS_obs_1[:, 2]
    k_2 = matter_PS_obs_2[:, 0]
    P_2 = matter_PS_obs_2[:, 1]
    h = 0.67
    err_2 = np.abs(matter_PS_obs_2[:, 2]-P_2)
    plt.plot(k_values, P_k_values/2/h, lw=2.5, label=r'$P_{\mathrm{M}}(k)$')
    plt.axvline(x=Hp_values[x_rad_is_mat]/MPC_to_meter/c/h, color='k', linestyle='--', lw=2.5, label=r'$k_{\mathrm{eq}}$')
    plt.errorbar(k_1, P_1, yerr=err_1, fmt='o', capsize=5, linewidth=0.7, markersize=4, label=r'SDSS DR7 LRG')
    plt.errorbar(k_2, P_2, yerr=err_2, fmt='o', capsize=5, linewidth=0.7, markersize=4, label=r'WMAP/ACT data')
    plt.xlim(1e-3, 0.4)
    plt.ylim(5e2, 5e4)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k\,\,[h/\mathrm{Mpc}]$')
    plt.ylabel(r'$P(k)\,\,\mathrm{[Mpc/}h]^3$')
    plt.legend()
    plt.title(r'Matter power spectrum $P(k)\,\,[\mathrm{Mpc}/h]^3$')
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/Mat_pow.pdf", format='pdf')



def C_ell(filename1, filename2):
    cosmic_variance = np.sqrt(2/(2*ell_values + 1))*c_ell_values
    plt.figure(figsize=(10, 6))
    plt.plot(ell_values, c_ell_values, lw=2.5, label=r'$\ell(\ell+1)C_\ell/(2\pi)$')
    x_l, y_l, lower_err_l, upper_err_l = read_data(filename1)
    x_h, y_h, lower_err_h, upper_err_h = read_data(filename2)
    plt.errorbar(x_l, y_l, yerr=[lower_err_l, upper_err_l], fmt='o', capsize=5, label=r'Low $\ell$ Planck data', linewidth=0.7, markersize=4, color = 'orange')
    plt.errorbar(x_h, y_h, yerr=[lower_err_h, upper_err_h], fmt='o', capsize=5, label=r'High $\ell$ Planck data', linewidth=0.7, markersize=4, color = 'red')
    plt.fill_between(ell_values, c_ell_values + cosmic_variance, c_ell_values - cosmic_variance, color = 'g', alpha = 0.3, label = 'Cosmic variance')
    plt.xscale('log')
    plt.xlim(1.95, 2000)
    plt.ylim(100, 6000)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell(\ell)/(2\pi)\,[\mu\mathrm{K}]^2$')
    plt.title(r'CMB power spectrum $\ell(\ell+1)C_\ell/(2\pi)$ as a function of $\ell$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/C_ell.pdf", format='pdf')

def C_ell_cheat(filename1, filename2):
    ell_values_cheat = ell_values**1.018
    c_ell_values_cheat = 1.13*c_ell_values*np.exp(-0.05*(ell_values/200)**(1.5))
    cosmic_variance_cheat = np.sqrt(2/(2*ell_values_cheat + 1))*c_ell_values_cheat
    plt.figure(figsize=(10, 6))
    plt.plot(ell_values_cheat, c_ell_values_cheat, lw=2.5, label=r'$\ell(\ell+1)C_\ell/(2\pi)$')
    x_l, y_l, lower_err_l, upper_err_l = read_data(filename1)
    x_h, y_h, lower_err_h, upper_err_h = read_data(filename2)
    plt.errorbar(x_l, y_l, yerr=[lower_err_l, upper_err_l], fmt='o', capsize=5, label=r'Low $\ell$ Planck data', linewidth=0.7, markersize=4, color = 'orange')
    plt.errorbar(x_h, y_h, yerr=[lower_err_h, upper_err_h], fmt='o', capsize=5, label=r'High $\ell$ Planck data', linewidth=0.7, markersize=4, color = 'red')
    plt.fill_between(ell_values_cheat, c_ell_values_cheat + cosmic_variance_cheat, c_ell_values_cheat - cosmic_variance_cheat, color = 'g', alpha = 0.3, label = 'Cosmic variance')
    plt.xscale('log')
    plt.xlim(1.95, 2000)
    plt.ylim(100, 6000)
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$\ell(\ell+1)C_\ell(\ell)/(2\pi)\,[\mu\mathrm{K}]^2$')
    plt.title(r'Adjusted CMB power spectrum $\ell(\ell+1)C_\ell/(2\pi)$ as a function of $\ell$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/C_ell_cheat.pdf", format='pdf')

def transfer_plot():
    plt.figure(figsize=(10, 6))
    theta_values_list = [theta_5_values, theta_20_values, theta_100_values, theta_200_values, theta_1000_values]
    ell_values = [5, 20, 100, 200, 1000]
    for theta_values, ell in zip(theta_values_list, ell_values):
        plt.plot(k_eta0_values, theta_values, label=rf'$\ell={ell}$', lw=1.2)
    plt.xlim(0, 800)
    plt.ylim(-0.015, 0.015)
    plt.xlabel(r'$k\eta_0$')
    plt.ylabel(r'$\Theta_\ell$')
    plt.title(r'$\Theta_\ell$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/transfer.pdf", format='pdf')

def integrand():
    plt.figure(figsize=(10, 6))
    theta_values_list = [theta_5_values, theta_20_values, theta_100_values, theta_200_values, theta_1000_values]
    ell_values = [5, 20, 100, 200, 1000]
    for theta_values, ell in zip(theta_values_list, ell_values):
        plt.plot(k_eta0_values, (theta_values**2 / k_eta0_values) * ell * (ell + 1), label=rf'$\ell={ell}$', lw=1.2)
    plt.xlim(1, 1500)
    plt.xlabel(r'$k\eta_0$')
    plt.xscale('log')
    plt.ylabel(r'$|\Theta_\ell(k\eta_0)|^2\ell(\ell+1)/(k\eta_0)$')
    plt.title(r'Integrand $|\Theta_\ell(k\eta_0)|^2/(k\eta_0)$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone4/Figures/integrand.pdf", format='pdf')

    


# Main function
def main():
    Mat_Pow()
    C_ell(planck_data_low, planck_data_high)
    planck_compare(planck_data_low)
    C_ell_cheat(planck_data_low, planck_data_high)
    transfer_plot()
    integrand()
    # plt.show()  # Uncomment to show plots, otherwise it just saves them to /Latex/Milestone4/Figures
main()
