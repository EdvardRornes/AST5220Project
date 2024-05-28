import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['axes.titlepad'] = 20 

font = {'family' : 'euclid',
        'weight' : 'bold',
        'size'   : 19}

matplotlib.rc('font', **font)

# Constants
MPC_to_meter   = 1/3.086e22

# I do not know why the following lines being the way they are is necessary but I could not run the script otherwise
# Get the current directory of the script
# current_directory    = os.path.dirname(os.path.abspath(__file__))

# Construct the paths to the files using the current directory
cosmology            = '../data/cosmology.txt'

perturbation_1       = '../data/perturbations_k0.1.txt'
perturbation_01      = '../data/perturbations_k0.01.txt'
perturbation_001     = '../data/perturbations_k0.001.txt'

# Extract from data
data_cosmology        = np.genfromtxt(cosmology)
data_perturbation_1   = np.genfromtxt(perturbation_1)
data_perturbation_01  = np.genfromtxt(perturbation_01)
data_perturbation_001 = np.genfromtxt(perturbation_001)

# Extract columns from cosmology
x_cos_values            = data_cosmology[:, 0]
eta_cos_values          = data_cosmology[:, 1]

# Extract columns from perturbation
x_values_1              = data_perturbation_1[:, 0]
T0_values_1             = data_perturbation_1[:, 1]
T1_values_1             = data_perturbation_1[:, 2]
T2_values_1             = data_perturbation_1[:, 3]
delta_cdm_values_1      = data_perturbation_1[:, 4]
delta_b_values_1        = data_perturbation_1[:, 5]
v_cdm_values_1          = data_perturbation_1[:, 6]
v_b_values_1            = data_perturbation_1[:, 7]
Phi_values_1            = data_perturbation_1[:, 8]
Psi_values_1            = data_perturbation_1[:, 9]
Source_T_values_1       = data_perturbation_1[:, 10]
Source_T_5_values_1     = data_perturbation_1[:, 11]
Source_T_50_values_1    = data_perturbation_1[:, 12]
Source_T_500_values_1   = data_perturbation_1[:, 13]
keta_values_1           = data_perturbation_1[:, 14]
Xe_values_1             = data_perturbation_1[:, 15]

x_values_01             = data_perturbation_01[:, 0]
T0_values_01            = data_perturbation_01[:, 1]
T1_values_01            = data_perturbation_01[:, 2]
T2_values_01            = data_perturbation_01[:, 3]
delta_cdm_values_01     = data_perturbation_01[:, 4]
delta_b_values_01       = data_perturbation_01[:, 5]
v_cdm_values_01         = data_perturbation_01[:, 6]
v_b_values_01           = data_perturbation_01[:, 7]
Phi_values_01           = data_perturbation_01[:, 8]
Psi_values_01           = data_perturbation_01[:, 9]
Source_T_values_01      = data_perturbation_01[:, 10]
Source_T_5_values_01    = data_perturbation_01[:, 11]
Source_T_50_values_01   = data_perturbation_01[:, 12]
Source_T_500_values_01  = data_perturbation_01[:, 13]
keta_values_01          = data_perturbation_01[:, 14]

x_values_001            = data_perturbation_001[:, 0]
T0_values_001           = data_perturbation_001[:, 1]
T1_values_001           = data_perturbation_001[:, 2]
T2_values_001           = data_perturbation_001[:, 3]
delta_cdm_values_001    = data_perturbation_001[:, 4]
delta_b_values_001      = data_perturbation_001[:, 5]
v_cdm_values_001        = data_perturbation_001[:, 6]
v_b_values_001          = data_perturbation_001[:, 7]
Phi_values_001          = data_perturbation_001[:, 8]
Psi_values_001          = data_perturbation_001[:, 9]
Source_T_values_001     = data_perturbation_001[:, 10]
Source_T_5_values_001   = data_perturbation_001[:, 11]
Source_T_50_values_001  = data_perturbation_001[:, 12]
Source_T_500_values_001 = data_perturbation_001[:, 13]
keta_values_001         = data_perturbation_001[:, 14]

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
                output = x_values_1[i]
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

hor_cross_1 = find_index(x_values_1, keta_values_1-1, -1, 0.01, 0, len(keta_values_1))
hor_cross_01 = find_index(x_values_01, keta_values_01-1, -1, 0.01, 0, len(keta_values_01))
hor_cross_001 = find_index(x_values_001, keta_values_001-1, -1, 0.01, 0, len(keta_values_001))
hor_cross_photon_1 = find_index(x_values_1, keta_values_1-1/np.sqrt(3), -1, 0.01, 0, len(keta_values_1))
hor_cross_photon_01 = find_index(x_values_01, keta_values_01-1/np.sqrt(3), -1, 0.01, 0, len(keta_values_01))
hor_cross_photon_001 = find_index(x_values_001, keta_values_001-1/np.sqrt(3), -1, 0.01, 0, len(keta_values_001))
recomb = idx_recombination = find_index(x_values_1, Xe_values_1, -1, 0.1, 0, len(x_values_1))

def data_test():
    plt.figure(figsize=(10, 6))
    plt.plot(x_cos_values, -1/3*np.cos(0.1*eta_cos_values/np.sqrt(3)*MPC_to_meter), lw=2.5, label=r'$\cos\left(\frac{k\eta}{\sqrt{3}}\right)$')
    plt.plot(x_values_1, T0_values_1 + Psi_values_1, lw=2.5, label=r'$\Theta_0+\Psi$')
    plt.xlim(-12, -7)    
    plt.xlabel(r'$x$')
    plt.title(r'Test of $\Theta_0+\Psi\propto\cos\left(\frac{k\eta}{\sqrt{3}}\right)$ with $k=0.1$ Mpc')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone3/Figures/test1.pdf", format='pdf')
    
def source():
    plt.figure(figsize=(10, 6))
    # plt.plot(x_values_1, Source_T_values_1, lw=2.5, label=r'$S(k=0.1)$')
    plt.plot(x_values_01, Source_T_values_01, lw=2.5, label=r'$S(k=0.01)$')
    # plt.plot(x_values_001, Source_T_values_001, lw=2.5, label=r'$S(k=0.001)$')
    plt.xlim(-8, 0)    
    plt.xlabel(r'$x$')
    plt.title(r'Test of $S(x)$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone3/Figures/source_function.pdf", format='pdf')

def gamma(k_value):
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))
    # I am very sorry that you had to witness this terrible code
    k_value_str = str(k_value)
    x_max = 0
    x_min = -13
    if (k_value == 0.1):
        y_min_delta = -3
        y_max_delta = 4
        y_min_v = -1.5
        y_max_v = 1.5
        use_x = x_values_1
        use_T0 = T0_values_1
        use_T1 = T1_values_1
        use_delta_b = delta_b_values_1
        use_v = v_b_values_1
        use_hor_cross = hor_cross_photon_1
    elif (k_value == 0.01):
        y_min_delta = 0.5
        y_max_delta = 4
        y_min_v = -0.75
        y_max_v = 0.75
        use_x = x_values_01
        use_T0 = T0_values_01
        use_T1 = T1_values_01
        use_delta_b = delta_b_values_01
        use_v = v_b_values_01
        use_hor_cross = hor_cross_photon_01
    elif (k_value == 0.001):
        y_min_delta = 0.75
        y_max_delta = 2.75
        y_min_v = -0.2
        y_max_v = 0.4
        use_x = x_values_001
        use_T0 = T0_values_001
        use_T1 = T1_values_001
        use_delta_b = delta_b_values_001
        use_v = v_b_values_001
        use_hor_cross = hor_cross_photon_001

    axs[0].plot(use_x, 4*use_T0, lw=2.5, label=r'$\delta_\gamma(k=' + k_value_str + r'/\mathrm{Mpc})$')
    axs[0].plot(use_x, use_delta_b, linestyle = '--', lw=2.5, label=r'$\delta_{\mathrm{B}}(k=' + k_value_str + r'/\mathrm{Mpc})$')
    axs[0].axvline(x=use_hor_cross, color='k', linestyle='-.', lw = 2.5, label='Horizon crossing')
    axs[0].axvline(x=recomb, linestyle='-.', lw = 2.5, label='Recombination')
    axs[0].set_xlim(x_min, x_max)
    axs[0].set_ylim(y_min_delta, y_max_delta)
    axs[0].set_xlabel(r'$x$')
    axs[0].set_ylabel(r'$\delta_\gamma$')
    axs[0].set_title(r'Density perturbation for photons $\delta_\gamma(x)=4\Theta_0(x)$')
    axs[0].legend()
    axs[0].grid(True)
    
    axs[1].plot(use_x, -3*use_T1, lw=2.5, label=r'$\delta_\gamma(k=' + k_value_str + r'/\mathrm{Mpc})$')
    axs[1].plot(use_x, use_v, linestyle = '--', lw=2.5, label=r'$v_b(k=' + k_value_str + r'/\mathrm{Mpc})$')
    axs[1].axvline(x=use_hor_cross, color='k', linestyle='-.', lw = 2.5, label='Horizon crossing')
    axs[1].axvline(x=recomb, linestyle='-.', lw = 2.5, label='Recombination')
    axs[1].set_xlim(x_min, x_max)
    axs[1].set_ylim(y_min_v, y_max_v)
    axs[1].set_xlabel(r'$x$')
    axs[1].set_ylabel(r'$v_\gamma$')
    axs[1].set_title(r'Velocity perturbation for photons $v_{\gamma}(x)=-3\Theta_1(x)$')
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    filename = 'gamma' + str(k_value)
    plt.savefig("../Latex/Milestone3/Figures/" + filename + ".pdf", format='pdf')

def gamma_all():
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    axs[0].plot(x_values_1, 4*T0_values_1, lw=2.5, label=r'$k=0.1/\mathrm{Mpc}$')
    axs[0].plot(x_values_01, 4*T0_values_01, lw=2.5, label=r'$k=0.01/\mathrm{Mpc}$')
    axs[0].plot(x_values_001, 4*T0_values_001, lw=2.5, label=r'$k=0.001/\mathrm{Mpc}$')
    axs[0].axvline(x=hor_cross_photon_1, color='b', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_photon_01, color='orange', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_photon_001, color='g', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[0].set_xlim(-13, 0)
    axs[0].set_ylim(-2.5, 4)
    axs[0].set_xlabel(r'$x$')
    axs[0].set_ylabel(r'$\delta_\gamma$')
    axs[0].set_title(r'Density perturbation for photons $\delta_\gamma(x)=4\Theta_0(x)$ at different modes $k$')
    axs[0].legend()
    axs[0].grid(True)
    
    axs[1].plot(x_values_1, -3*T1_values_1, lw=2.5, label=r'$k=0.1/\mathrm{Mpc}$')
    axs[1].plot(x_values_01, -3*T1_values_01, lw=2.5, label=r'$k=0.01/\mathrm{Mpc}$')
    axs[1].plot(x_values_001, -3*T1_values_001, lw=2.5, label=r'$k=0.001/\mathrm{Mpc}$')
    axs[1].axvline(x=hor_cross_photon_1, color='b', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_photon_01, color='orange', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_photon_001, color='g', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[1].set_xlim(-13, 0)
    axs[1].set_ylim(-1.5, 1.5)
    axs[1].set_xlabel(r'$x$')
    axs[1].set_ylabel(r'$v_\gamma$')
    axs[1].set_title(r'Velocity perturbation for photons $v_{\gamma}(x)=-3\Theta_1(x)$ at different modes $k$')
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    plt.savefig("../Latex/Milestone3/Figures/gamma.pdf", format='pdf')

    

def b_cdm():
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    axs[0].plot(x_values_1, np.abs(delta_b_values_1), lw=2.5, label=r'$|\delta_{\mathrm{B}}|(0.1)$')
    axs[0].plot(x_values_01, delta_b_values_01, lw=2.5, label=r'$\delta_{\mathrm{B}}(0.01)$')
    axs[0].plot(x_values_001, delta_b_values_001, lw=2.5, label=r'$\delta_{\mathrm{B}}(0.001)$')
    axs[0].plot(x_values_1, delta_cdm_values_1, linestyle='--', lw=2.5, label=r'$\delta_{\mathrm{CDM}}(0.1)$')
    axs[0].plot(x_values_01, delta_cdm_values_01, linestyle='--', lw=2.5, label=r'$\delta_{\mathrm{CDM}}(0.01)$')
    axs[0].plot(x_values_001, delta_cdm_values_001, linestyle='--', lw=2.5, label=r'$\delta_{\mathrm{CDM}}(0.001)$')
    axs[0].axvline(x=hor_cross_1, color = 'b', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_01, color = 'orange', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_001, color = 'g', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[0].set_xlim(-13, 0)
    axs[0].set_ylim(1e-1, 1e5)
    axs[0].set_yscale('log')
    axs[0].set_xlabel(r'$x$')
    axs[0].set_ylabel(r'$\delta_i$')
    axs[0].set_title(r'Density perturbation for baryons $\delta_{\mathrm{B}}(k,x)$ and CDM $\delta_{\mathrm{CDM}}(k,x)$')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(x_values_1, np.abs(v_b_values_1), lw=2.5, label=r'$|v_{\mathrm{B}}|(0.1)$')
    axs[1].plot(x_values_01, v_b_values_01, lw=2.5, label=r'$v_{\mathrm{B}}(0.01)$')
    axs[1].plot(x_values_001, v_b_values_001, lw=2.5, label=r'$v_{\mathrm{B}}(0.001)$')
    axs[1].plot(x_values_1, v_cdm_values_1, linestyle='--', lw = 2.5, label=r'$v_{\mathrm{CDM}}(0.1)$')
    axs[1].plot(x_values_01, v_cdm_values_01, linestyle='--', lw = 2.5, label=r'$v_{\mathrm{CDM}}(0.01)$')
    axs[1].plot(x_values_001, v_cdm_values_001, linestyle='--', lw = 2.5, label=r'$v_{\mathrm{CDM}}(0.001)$')
    axs[1].axvline(x=hor_cross_1, color = 'b', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_01, color = 'orange', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_001, color = 'g', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[1].set_xlim(-13, 0)
    axs[1].set_ylim(1e-3, 1e2)
    axs[1].set_yscale('log')
    axs[1].set_xlabel(r'$x$')
    axs[1].set_ylabel(r'$v_i$')
    axs[1].set_title(r'Velocity perturbation for baryons $v_{\mathrm{B}}(k,x)$ and CDM $v_{\mathrm{CDM}}(k,x)$')
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    plt.savefig("../Latex/Milestone3/Figures/baryon_CDM.pdf", format='pdf')


def Phi_Psi():
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    axs[0].plot(x_values_1, Phi_values_1, lw=2.5, label=r'$k=0.1$')
    axs[0].plot(x_values_01, Phi_values_01, lw=2.5, label=r'$k=0.01$')
    axs[0].plot(x_values_001, Phi_values_001, lw=2.5, label=r'$k=0.001$')
    axs[0].axvline(x=hor_cross_1, color = 'b', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_01, color = 'orange', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=hor_cross_001, color = 'g', linestyle='-.', lw = 2.5)
    axs[0].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[0].set_xlim(-15, 0)
    axs[0].set_xlabel(r'$x$')
    axs[0].set_ylabel(r'$\Phi$')
    axs[0].set_title(r'Gravitational potential $\Phi(x)$ at different modes $k$')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(x_values_1, Phi_values_1 + Psi_values_1, lw=2.5, label=r'$k=0.1$')
    axs[1].plot(x_values_01, Phi_values_01 + Psi_values_01, lw=2.5, label=r'$k=0.01$')
    axs[1].plot(x_values_001, Phi_values_001 + Psi_values_001, lw=2.5, label=r'$k=0.001$')
    axs[1].axvline(x=hor_cross_1, color = 'b', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_01, color = 'orange', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=hor_cross_001, color = 'g', linestyle='-.', lw = 2.5)
    axs[1].axvline(x=recomb, color='k', linestyle='-.', lw = 2.5, label='Recombination')
    axs[1].set_xlim(-15, 0)
    axs[1].set_xlabel(r'$x$')
    axs[1].set_ylabel(r'$\Phi+\Psi$')
    axs[1].set_title(r'Sum of the gravitational potentials $\Phi+\Psi$ at different modes $k$')
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    plt.savefig("../Latex/Milestone3/Figures/Phi_Psi.pdf", format='pdf')

def T2():
    plt.figure(figsize=(10, 6))
    plt.plot(x_values_1, T2_values_1, lw=2.5, label=r'$k=0.1$')
    plt.plot(x_values_01, T2_values_01, lw=2.5, label=r'$k=0.01$')
    plt.plot(x_values_001, T2_values_001, lw=2.5, label=r'$k=0.001$')
    plt.axvline(x=recomb, color='k', linestyle='--', lw = 2.5, label='Recombination')
    plt.axvline(x=hor_cross_1, color = 'b', linestyle='-.', lw = 2.5)
    plt.axvline(x=hor_cross_01, color = 'orange', linestyle='-.', lw = 2.5)
    plt.axvline(x=hor_cross_001, color = 'g', linestyle='-.', lw = 2.5)
    plt.xlim(-13, 0)    
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\Theta_2$')
    plt.title(r'Photon multiple $\Theta_2(x)$ at different modes $k$')
    plt.legend()
    plt.grid(True)
    plt.savefig("../Latex/Milestone3/Figures/Theta2.pdf", format='pdf')


# Main function
def main():
    source()
    data_test()
    gamma(0.1)
    gamma(0.01)
    gamma(0.001)
    gamma_all()
    b_cdm()
    Phi_Psi()
    T2()
    # plt.show()  # Uncomment to show plots, otherwise it just saves them to /Latex/Milestone3/Figures


main()
