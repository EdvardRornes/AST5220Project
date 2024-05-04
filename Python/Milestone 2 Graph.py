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
Hkm_per_parsec = 10**5*MPC_to_meter
Gigayear       = 10**9*365*24*3600 # Gigayear to seconds
H0             = 0.67*10**5*MPC_to_meter

# I do not know why the following lines being the way they are is necessary but I could not run the script otherwise
# Get the current directory of the script
# current_directory    = os.path.dirname(os.path.abspath(__file__))

# Construct the paths to the files using the current directory
recombination      = '../data/recombination.txt'

# Extract from data
data_recombination = np.genfromtxt(recombination)

# Extract columns from cosmology
x_values                 = data_recombination[:, 0]
Xe_values                = data_recombination[:, 1]
Xe_saha_values           = data_recombination[:, 2]
ne_values                = data_recombination[:, 3]
tau_values               = data_recombination[:, 4]
dtaudx_values            = data_recombination[:, 5]
ddtauddx_values          = data_recombination[:, 6]
g_tilde_values           = data_recombination[:, 7]
max_g_tilde              = np.max(data_recombination[:, 7])
dgdx_norm_tilde_values   = data_recombination[:, 8]/np.max(data_recombination[:, 8])*max_g_tilde
ddgddx_norm_tilde_values = data_recombination[:, 9]/np.max(abs(data_recombination[:, 9]))*max_g_tilde
sound_horizon_values     = data_recombination[:, 10]
t_values                 = data_recombination[:, 11]
z_values                 = data_recombination[:, 12]



# Index finder function
def find_index(array, type, threshold, start, end):
    for i in range(start, end):
        if type == 0:
            if array[i] == threshold:
                output = x_values[i]
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
                output = x_values[i]
                return output
            else:
                output = None
        elif type == -1:
            if np.abs(array[i]) <= threshold:
                output = x_values[i]
                return output
        elif type == -2:
            if np.abs(array[i]) <= threshold:
                output = i
                return output
            else:
                output = None
        else:
            print('Not valid condition type')
            break


def tau(x_last_scattering):
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, tau_values, lw=2.5, label=r'$\tau(x)$')
    plt.plot(x_values, -dtaudx_values, lw=2.5, label=r'$-\tau^{\prime}(x)$')
    plt.plot(x_values, ddtauddx_values, lw=2.5, label=r'$\tau^{\prime\prime}(x)$')
    plt.axvline(x=x_last_scattering, color='k', linestyle='--', lw = 1.5, label='Last scattering')
    plt.xlim(-10, 0)
    plt.ylim(1e-7, 1e5)
    plt.xlabel(r'$x$')
    plt.ylabel(r'Optical depth')
    plt.yscale('log')
    plt.title(r'Optical depth $\tau(x)$ and its derivatives')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone2/Figures/tau_and_derivs.pdf", format='pdf')

def ne():
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, ne_values, lw=2.5, label=r'$n_e(x)$')
    plt.xlim(-10, 0)
    plt.ylim(1e-6, 1e14)
    plt.xlabel(r'$x$')
    plt.ylabel(r'Number density')
    plt.yscale('log')
    plt.title(r'Electron number density $n_e(x)$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone2/Figures/electron_number_density.pdf", format='pdf')

def Xe(x_recombination):
    # Find freeze out Xe value, i.e. Xe(x=0)
    idx_freeze_out = find_index(x_values, -2, 1e-5, 0, len(x_values))
    Xe_freeze_out = Xe_values[idx_freeze_out]
    print(Xe_freeze_out)
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, Xe_values, lw=2.5, label=r'$X_e(x)$')
    plt.plot(x_values, Xe_saha_values, lw=2.5, label=r'$X_e^\mathrm{Saha}(x)$')
    plt.axvline(x=x_recombination, color='k', linestyle='--', lw = 1.5, label='Recombination')
    plt.axhline(y=Xe_freeze_out, color='r', linestyle='--', lw = 1.5, label='Freeze out')
    plt.xlim(-7.8, -4)
    plt.ylim(1e-4, 1.5)
    plt.xlabel(r'$x$')
    plt.yscale('log')
    plt.title(r'Fraction of free electrons $X_e(x)$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone2/Figures/Xe.pdf", format='pdf')

def g_tilde(x_last_scattering):
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, g_tilde_values, lw=2.5, label=r'$\tilde{g}(x)$')
    plt.plot(x_values, dgdx_norm_tilde_values, lw=2.5, label=r'$\tilde{g}^{\prime}(x)$')
    plt.plot(x_values, ddgddx_norm_tilde_values, lw=2.5, label=r'$\tilde{g}^{\prime\prime}(x)$')
    plt.axvline(x=x_last_scattering, color='k', linestyle='--', lw = 1.5, label='Last scattering')
    plt.xlim(-7.25, -6.5)
    plt.ylim(-1.1*max_g_tilde, 1.1*max_g_tilde)
    plt.xlabel(r'$x$')
    plt.title(r'Visibility function $\tilde{g}(x)$ its derivatives normalized')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone2/Figures/gtilde.pdf", format='pdf')

    # Make a table for the values of the various time parameters in the various relevant cosmic events
def make_table():
    # Find the values where recombination ocurred
    idx_recombination = find_index(Xe_values, -2, 0.1, 0, len(x_values))
    print(idx_recombination)
    x_recombination_str = str(round(x_values[idx_recombination], 6))
    z_recombination_str = str(round(z_values[idx_recombination]))
    t_recombination_str = str(round(t_values[idx_recombination]/Gigayear*10**9))
    s_recombination_str = str(round(sound_horizon_values[idx_recombination]*MPC_to_meter, 1))

    # Find the values where last scattering/decoupling ocurred
    idx_last_scattering = find_index(g_tilde_values, 2, np.max(g_tilde_values), 0, len(x_values))
    print(idx_last_scattering)
    x_last_scattering_str = str(round(x_values[idx_last_scattering], 6))
    z_last_scattering_str = str(round(z_values[idx_last_scattering]))
    t_last_scattering_str = str(round(t_values[idx_last_scattering]/Gigayear*10**9))
    s_last_scattering_str = str(round(sound_horizon_values[idx_last_scattering]*MPC_to_meter, 1))

    table_data = [
        ['x'         , x_recombination_str, x_last_scattering_str],
        ['z(x)'      , z_recombination_str, z_last_scattering_str],
        ['t(x) [yr]', t_recombination_str, t_last_scattering_str],
        ['s(x) [Mpc]', s_recombination_str, s_last_scattering_str]
    ]


    # Create the table
    table = tabulate(table_data, headers=["Event", "Recombination", "Last scattering/decoupling"], tablefmt="simple_grid")

    print(table)

    return x_values[idx_recombination], x_values[idx_last_scattering]
      

# Main function
def main():
    x_recombination, x_last_scattering = make_table()
    Xe(x_recombination)
    tau(x_last_scattering)
    ne()
    g_tilde(x_last_scattering)
    # plt.show()

main()
