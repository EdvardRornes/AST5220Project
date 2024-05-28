# Cosmology II Project

This repository is built on a C++ templates from [website](https://github.com/HAWinther/AST5220-Cosmology) for making an Einstein-Boltzmann solver (a CAMB like code). This is the project for the course AST5220 "Cosmology II" at ITA University of Oslo. 

### Website
All relevant information background for this project can also be found at [website](https://cmb.wintherscoming.no/).

## Structure

### Background Cosmology

In this section we consider the evolution of the uniform background of the universe in the Lambda-Cold-Dark-Matter ($\Lambda$CDM) model which is considered to be the present day standard model of cosmology. The main goal of this section is to solve the unperturbed background evolution of the universe by using known cosmological parameters and comparing this to observational supernova data.

### Recombination History

Next we look at the recombination history of the universe. This is when baryons, mainly protons and electrons, went from being ionized to forming neutral atoms once the energy of the photons dropped below $13.6\,$eV. As a result, photons during this time decoupled from the thermal equilibrium of the universe and are what we now detect as the CMB photons. As this event is tightly related to the free mean path of photons, the goal of this section is to compute the optical depth $\tau$, the visibility function $\tilde g$ and their derivatives as these will be vital for the next sections. In this section we only consider the formation of hydrogen and neglect the existence of any heavier atoms.

### Perturbations

The topic of this section is to see how small quantum fluctuations from the inflationary period caused small perturbations to the baryon, photon and dark matter fluid in the early which then grew into large structures formations we see today. Since we have determined various quantities of the background cosmology in the previous sections, all we need to do is to perturb the background and determine suitable initial conditions. To do this we write the perturbed flat FRLW metric in the Newtonian gauge. Note that from this section and onwards we set $N_\text{eff}=0$, completely neglecting neutrinos from hereon.

### Power Spectrum

In this section we compute the CMB power spectrum using the line-of-sight integration technique and compare the theoretical prediction with observables.

## Prerequisites

The code uses GSL, see [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac see [this](https://cmb.wintherscoming.no/about.php#GSL) for how to install it.

## Running the C++ code

Build:

    make cmb

Run:

 ./cmb

Clean:

    make clean

## Running the Python code

All plots are saved in /Latex/MilestoneX/Figures. To see the graphs, uncomment the plt.show() for the desired section at the bottom of /Python/MilestoneX_1_Graph.py
