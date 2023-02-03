import numpy as np
import subprocess as sp


def get_MO_MATRIX():
    
    N_COLUMNS = 5 # I think always 5 columns in FCHK files

    N_ORBITAL_ENERGIES = int(sp.check_output(\
        "grep 'Alpha Orbital Energies' geometry.fchk | awk '{print $6}'", shell=True))
    
    N_ROWS = int(np.ceil(N_ORBITAL_ENERGIES / N_COLUMNS))
    sp.call(f"grep 'Alpha Orbital Energies' geometry.fchk -A {N_ROWS} | tail -n {N_ROWS} > MO_ENERGY.dat", shell=True)
    MO_ENERGY = np.array([ float(it) for line in open("MO_ENERGY.dat","r").readlines() for it in line.split() ])

    N_MO_COEFFS = N_ORBITAL_ENERGIES**2
    N_ROWS = int(np.ceil(N_MO_COEFFS / N_COLUMNS))
    sp.call(f"grep 'Alpha MO coefficients' geometry.fchk -A {N_ROWS} | tail -n {N_ROWS} > MO_COEFFS.dat", shell=True)
    MO_COEFFS = np.array([ float(it) for line in open("MO_COEFFS.dat","r").readlines() for it in line.split() ])

    MO_COEFFS = MO_COEFFS.reshape(( N_ORBITAL_ENERGIES, N_ORBITAL_ENERGIES ))

    return MO_ENERGY, MO_COEFFS

if ( __name__ == "__main__" ):
    MO_ENERGY, MO_COEFFS = get_MO_MATRIX()