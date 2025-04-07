import numpy as np
from matplotlib import pyplot as plt

## Trajectories
NTRAJ = 200

## For Plotting
SIGMA = 0.0001 # eV
EMIN  = 4.0 # eV
EMAX  = 7.0 # eV
NPTS  = 1000

## Initialize Variables
TIME            = np.loadtxt("TRAJ/traj-0/MD_OUTPUT/PES.dat")[:,0]
NSTEPS, NSTATES = np.loadtxt("TRAJ/traj-0/MD_OUTPUT/PES.dat")[:,1:].shape # Skip time column
print(f"Number of Steps: {NSTEPS}")
print(f"Number of States: {NSTATES}")

## Read the data
ENERGY = np.zeros((NTRAJ, NSTEPS, NSTATES))
DIPOLE = np.zeros((NTRAJ, NSTEPS, NSTATES, 3))
for traj in range(NTRAJ):
    try:
        print(f"Reading traj = {traj} of {NTRAJ}")
        ENERGY[traj,:,:]   = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/PES.dat")[:NSTEPS,1:] # Skip time column
        DIPOLE[traj,:,:,0] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_X.dat")[:NSTEPS,1:] # Skip time column
        DIPOLE[traj,:,:,1] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_Y.dat")[:NSTEPS,1:] # Skip time column
        DIPOLE[traj,:,:,2] = np.loadtxt(f"TRAJ/traj-{traj}/MD_OUTPUT/S0_Sn_Dipoles_Z.dat")[:NSTEPS,1:] # Skip time column
    except ValueError:
        print(f"Trajectory {traj} is not complete. Skipping...")

EGRID = np.linspace(EMIN, EMAX, NPTS)
ABS_G = np.zeros(NPTS)
ABS_L = np.zeros(NPTS)
for traj in range(NTRAJ):
    print(f"Working on traj = {traj} of {NTRAJ}.")
    for step in range(NSTEPS):
        E0n = ENERGY[traj,step,1:] - ENERGY[traj,step,0]
        dE  = EGRID[:,None] - E0n[None,:]
        OSC = (2/3) * E0n * np.einsum("jd,jd->j", DIPOLE[traj,step,1:,:], DIPOLE[traj,step,1:,:])
        ABS_G[:] += np.sum( OSC[None,:] * np.exp(-dE[:,:]**2/2/SIGMA**2), axis=1 )
        ABS_L[:] += np.sum( OSC[None,:] * SIGMA**2 / (SIGMA**2 + dE[:,:]**2), axis=1 )

ABS_G /= np.max(ABS_G) # NSTEPS * NTRAJ
ABS_L /= np.max(ABS_L) # NSTEPS * NTRAJ

plt.plot(EGRID, ABS_G, "-", c="black", lw=3, label="Gaussian")
plt.plot(EGRID, ABS_L, "--", c="red", lw=2,  label="Lorentzian")
plt.legend()
plt.xlim(EMIN,EMAX)
plt.ylim(0)
plt.xlabel("Energy (eV)", fontsize=15)
plt.ylabel("Absorption (Arb. Units)", fontsize=15)
plt.savefig("Absorption_FreqDomain.jpg", dpi=300)


