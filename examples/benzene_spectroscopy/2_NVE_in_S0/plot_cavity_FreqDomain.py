import numpy as np
from matplotlib import pyplot as plt
import sys

## Trajectories
NTRAJ = 50 # Number of trajectories

## Cavity
lam_collective = 0.02 # Collective cavity coupling (a.u.) -- wc A0 = sqrt(wc/2) lam_c = gc / mu
wc            = 5.425 / 27.2114 # Cavity frequency (a.u.)
if ( len(sys.argv) >= 2 ):
    NMOL          = int( sys.argv[1] ) # Number of molecules collectively coupled
NEL           = 10 # Number of electronic states per molecule
NPOL          = 1 + NMOL*(NEL-1) + 1 # GS, MATTER (mol+el), cavity
lam_single    = lam_collective / np.sqrt(NMOL) # Single molecule coupling (a.u.)
pol           = np.array([1,1,1]) # Polarization vector (a.u.)
pol           = pol / np.linalg.norm(pol)

NAVERAGES = NTRAJ - NMOL + 1 # Sliding averages

## For Plotting
SIGMA = 0.0001 # eV
EMIN  = 4.0 # eV
EMAX  = 7.0 # eV
NPTS  = 1000
EGRID = np.linspace(EMIN, EMAX, NPTS)

## Initialize Variables
TIME            = np.loadtxt("TRAJ/traj-0/MD_OUTPUT/PES.dat")[:,0]
NSTEPS, NSTATES = np.loadtxt("TRAJ/traj-0/MD_OUTPUT/PES.dat")[:,1:].shape # Skip time column
print(f"Number of Steps: {NSTEPS}")
print(f"Number of States: {NSTATES}")

## Read the data
E_MOL = np.zeros((NTRAJ, NSTEPS, NSTATES))
DIP_MOL = np.zeros((NTRAJ, NSTEPS, NSTATES, 3))
for traj in range(NTRAJ):
    try:
        print(f"Reading traj = {traj} of {NTRAJ}")
        E_MOL[traj,:,:]     = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/PES.dat")[:NSTEPS,1:] / 27.2114 # Skip time column
        DIP_MOL[traj,:,:,0] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_X.dat")[:NSTEPS,1:] # Skip time column
        DIP_MOL[traj,:,:,1] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_Y.dat")[:NSTEPS,1:] # Skip time column
        DIP_MOL[traj,:,:,2] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_Z.dat")[:NSTEPS,1:] # Skip time column
    except ValueError:
        print(f"Trajectory {traj} is not complete. Skipping...")

DIP_MOL = DIP_MOL.dot(pol) # Project onto polarization vector
E_JC = np.zeros((NAVERAGES, NSTEPS, NPOL))
PHOT = np.zeros((NAVERAGES, NSTEPS, NPOL))
PDIP = np.zeros((NAVERAGES, NSTEPS, NPOL))
for avei in range(NAVERAGES):
    print(f"Calculating Hamiltonian {avei} of {NAVERAGES} averages with {NTRAJ} trajectories")
    for ti in range( NSTEPS ):
        H_JC = np.zeros((NPOL, NPOL))
        # Diagonal Energies
        H_JC[0,0] = np.sum( E_MOL[:NMOL+avei,ti,0] )
        for A in range( NMOL ):
            for j in range( NEL-1 ):
                H_JC[A*(NEL-1)+1+j,A*(NEL-1)+1+j] = H_JC[0,0] + E_MOL[A+avei,ti,j+1] - E_MOL[A+avei,ti,0]
        H_JC[-1,-1] = H_JC[0,0] + wc

        # LM Coupling
        for A in range( NMOL ):
            for j in range( NEL-1 ):
                H_JC[A*(NEL-1)+1+j,-1] = np.sqrt(wc/2) * lam_single * DIP_MOL[A+avei,ti,j+1]
                H_JC[-1,A*(NEL-1)+1+j] = H_JC[A*(NEL-1)+1+j,-1]
        e,u = np.linalg.eigh(H_JC)
        E_JC[avei,ti,:] = e * 27.2114 # Convert to eV
        PHOT[avei,ti,:] = u[-1,:].conj() * u[-1,:] # Photon Population
        PDIP_OP    = np.zeros((NPOL,NPOL))
        for A in range( NMOL ):
            for j in range( NEL-1 ):
                PDIP_OP[0,A*(NEL-1)+1+j] = DIP_MOL[A+avei,ti,j+1] # < G, 0 | \hat{\mu} | E, 0 >
                PDIP_OP[A*(NEL-1)+1+j,0] = DIP_MOL[A+avei,ti,j+1] # < G, 0 | \hat{\mu} | E, 0 >
        #PDIP[avei,ti,:] = np.sum( u.conj() @ PDIP_OP @ u, axis=1 ) # Dipole Population
        PDIP[avei,ti,:] = np.einsum("xj,xy,yk->jk", u.conj(), PDIP_OP, u )[0,:] # Dipole Population
    sys.stdout.flush()

E_JC = E_JC.reshape( NAVERAGES*NSTEPS, NPOL )
PHOT = PHOT.reshape( NAVERAGES*NSTEPS, NPOL )
PDIP = PDIP.reshape( NAVERAGES*NSTEPS, NPOL )

ABS_G = np.zeros(NPTS)
ABS_L = np.zeros(NPTS)
TM_G  = np.zeros(NPTS)
TM_L  = np.zeros(NPTS)
for avei in range(E_JC.shape[0]):
    if ( avei%1000 == 0 ):
        print(f"Calculating spectra for {avei} of {E_JC.shape[0]} averages")
    E0n = E_JC[avei,1:] - E_JC[avei,0] # Ground-to-excited E_MOL difference
    dE  = EGRID[:,None] - E0n[None,:] # Shifted Gaussian/Lorentzian location
    GFUNC = np.exp(-dE[:,:]**2/2/SIGMA**2)
    LFUNC = SIGMA**2 / (SIGMA**2 + dE[:,:]**2)
    ABS_G[:] += np.sum( 2 * E0n[None,:] * PDIP[avei,1:][None,:]**2 * GFUNC, axis=1 ) # A Gaussian
    ABS_L[:] += np.sum( 2 * E0n[None,:] * PDIP[avei,1:][None,:]**2 * LFUNC, axis=1 ) # A Lorentzian
    TM_G[:]  += np.sum( PHOT[avei,1:][None,:] * GFUNC, axis=1 ) # A Gaussian
    TM_L[:]  += np.sum( PHOT[avei,1:][None,:] * LFUNC, axis=1 ) # A Lorentzian
    sys.stdout.flush()

ABS_G /= np.max(ABS_G) # NSTEPS * NTRAJ
ABS_L /= np.max(ABS_L) # NSTEPS * NTRAJ
TM_G /= np.max(TM_G) # NSTEPS * NTRAJ
TM_L /= np.max(TM_L) # NSTEPS * NTRAJ

plt.plot(EGRID, ABS_G, "-", c="black", lw=3, label="Gaussian")
plt.plot(EGRID, ABS_L, "--", c="red", lw=2,  label="Lorentzian")
plt.legend()
plt.xlim(EMIN,EMAX)
plt.ylim(0)
plt.xlabel("E_MOL (eV)", fontsize=15)
plt.ylabel("Absorption (Arb. Units)", fontsize=15)
plt.savefig("ABS_Cavity_FreqDomain_WC_%1.5f_LAM_%1.5f_NMOL_%d.jpg" % (wc*27.2114, lam_collective, NMOL), dpi=300)
plt.clf()
plt.close()

plt.plot(EGRID, TM_G, "-", c="black", lw=3, label="Gaussian")
plt.plot(EGRID, TM_L, "--", c="red", lw=2,  label="Lorentzian")
plt.legend()
plt.xlim(EMIN,EMAX)
plt.ylim(0)
plt.xlabel("E_MOL (eV)", fontsize=15)
plt.ylabel("Transmission (Arb. Units)", fontsize=15)
plt.savefig("TM_Cavity_FreqDomain_WC_%1.5f_LAM_%1.5f_NMOL_%d.jpg" % (wc*27.2114, lam_collective, NMOL), dpi=300)
plt.clf()
plt.close()

