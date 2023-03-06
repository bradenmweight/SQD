# Documentation for SQD Package

For questions, concerns, or bugs, please send an email to Braden M. Weight (<bweight@ur.rochester.edu>). The relevant parmeters can be found on the [parameters][PARAMETERS] page.

[PARAMETERS]: <https://bradenmweight.github.io/SQD/read.html?filename=Parameters.md>

# Table of Contents
- [Documentation for SQD Package](#documentation-for-sqd-package)
- [Table of Contents](#table-of-contents)
  - [**Mixed Quantum-Classical and Semi-classical Dynamics Methods**](#mixed-quantum-classical-and-semi-classical-dynamics-methods)
    - [**Born-Oppenheimer Molecular Dyanmics**](#born-oppenheimer-molecular-dyanmics)
    - [**Ehrenfest**](#ehrenfest)
    - [**Spin-LSC**](#spin-lsc)
    - [**Global Flux Surface Hopping**](#global-flux-surface-hopping)
    - [**Classical Path Approximation**](#classical-path-approximation)
  - [**Electronic Time-Propagation**](#electronic-time-propagation)
    - [**Quasi-diabatic (Local Diabatization) Propagation**](#quasi-diabatic-local-diabatization-propagation)
    - [**Second-order symplectic ("Velocity-Verlet-like")**](#second-order-symplectic-velocity-verlet-like)
    - [**Explicit Fourth-order Runge-Kutta**](#explicit-fourth-order-runge-kutta)
  - [**Statistical Ensembles**](#statistical-ensembles)
    - [**Canonical Ensemble (NVT)**](#canonical-ensemble-nvt)
    - [**Micro-canonical Ensemble (NVE)**](#micro-canonical-ensemble-nve)
  - [**Initialization of Nuclear Velocity**](#initialization-of-nuclear-velocity)
    - [**Read from File**](#read-from-file)
    - [**Draw from Maxwell-Boltzmann Distribution**](#draw-from-maxwell-boltzmann-distribution)
    - [**Set Velocities to Zero**](#set-velocities-to-zero)
  - [**Suggested Workflow**](#suggested-workflow)
    - [**Ground State Dynamics**](#ground-state-dynamics)
    - [**Sample Nuclear Distribution**](#sample-nuclear-distribution)
    - [**Initialize Non-adiabatic Trajectories**](#initialize-non-adiabatic-trajectories)
  - [**References**](#references)

## **Mixed Quantum-Classical and Semi-classical Dynamics Methods**
### **Born-Oppenheimer Molecular Dyanmics**
${\bf F} (t) = \langle S_{INIT} | \boldsymbol{\nabla} \hat{H}_\mathrm{el} | S_{INIT} \rangle $
### **Ehrenfest**
### **Spin-LSC**
### **Global Flux Surface Hopping**
### **Classical Path Approximation**
$${\bf F} (t) = \langle S_{0} | \boldsymbol{\nabla} \hat{H}_\mathrm{el} | S_{0} \rangle $$
## **Electronic Time-Propagation**
### **Quasi-diabatic (Local Diabatization) Propagation**
Math looks like this: 
$$S^{\dag} \times z_{\mu} \rightarrow z_{\mu}$$ 
### **Second-order symplectic ("Velocity-Verlet-like")**
### **Explicit Fourth-order Runge-Kutta**
## **Statistical Ensembles**
### **Canonical Ensemble (NVT)**
### **Micro-canonical Ensemble (NVE)**
## **Initialization of Nuclear Velocity**
### **Read from File**
### **Draw from Maxwell-Boltzmann Distribution**
### **Set Velocities to Zero**
## **Suggested Workflow**
### **Ground State Dynamics**
### **Sample Nuclear Distribution**
### **Initialize Non-adiabatic Trajectories**
## **References**
1. Mark Tuckerman, [Statistical Mechanics: Theory and Molecular Simulation](https://books.google.com/books?id=Lo3Jqc0pgrcC)
2. Mandal, Yamijala, and Huo, [Quasi-Diabatic Representation for Nonadiabatic Dynamics Propagation](https://pubs.acs.org/doi/10.1021/acs.jctc.7b01178), *J. Chem. Theory Comput.* 2018, 14, 4, 1828–1840 
3. Mandal, Sandoval, Shakib, and Huo, [Quasi-Diabatic Propagation Scheme for Direct Simulation of Proton-Coupled Electron Transfer Reaction](https://pubs.acs.org/doi/10.1021/acs.jpca.9b00077), *J. Phys. Chem. A* 2019, 123, 12, 2470–2482
4. Sandoval, Mandal, and Huo, [Symmetric quasi-classical dynamics with quasi-diabatic propagation scheme](https://aip.scitation.org/doi/full/10.1063/1.5036787), *J. Chem. Phys.* 149, 044115 (2018)
5. Weight, Mandal, and Huo, [Ab initio symmetric quasi-classical approach to investigate molecular Tully models](https://aip.scitation.org/doi/10.1063/5.0061934), *J. Chem. Phys.* 155, 084106 (2021)
6. Runeson and Richardson, [Spin-mapping approach for nonadiabatic molecular dynamics](https://aip.scitation.org/doi/10.1063/1.5100506), *J. Chem. Phys.* 151, 044119 (2019)
7. Runeson and Richardson, [Generalized spin mapping for quantum-classical dynamics](https://aip.scitation.org/doi/full/10.1063/1.5143412), *J. Chem. Phys.* 152, 084110 (2020)
8. Mannouch and Richardson, [A partially linearized spin-mapping approach for nonadiabatic dynamics. I. Derivation of the theory](https://aip.scitation.org/doi/full/10.1063/5.0031168), *J. Chem. Phys.* 153, 194109 (2020)
9.  Mannouch and Richardson, [A partially linearized spin-mapping approach for nonadiabatic dynamics. II. Analysis and comparison with related approaches](https://aip.scitation.org/doi/full/10.1063/5.0031173), *J. Chem. Phys.* 153, 194110 (2020)