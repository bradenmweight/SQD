import numpy as np
from matplotlib import pyplot as plt

TEMP_LIST = np.array([10,300,1000])
NTRAJ     = 1000 # 3000
NSTATES   = 1

NMOL      = 100
lambda_collective = 0.02
lambda_single     = lambda_collective / np.sqrt( NMOL )
NPOL      = NMOL + 1

NAVERAGES = NTRAJ - NMOL + 1
iSTATE    = 1
pol       = np.array([1,1,1])
pol       = pol / np.linalg.norm(pol)

ENERGY = np.zeros( (len(TEMP_LIST), NTRAJ) )
DIPOLE = np.zeros( (len(TEMP_LIST), NTRAJ) )
for Ti,T in enumerate( TEMP_LIST ):
    print(Ti, T)
    for traj in range( NTRAJ ):
        e = np.loadtxt(f"{T}K/2_TDDFT/TRAJ/traj-{traj}/PLOTS_DATA/ADIABATIC_ENERGIES_RPA.dat")
        ENERGY[Ti,traj] = e[iSTATE] - e[0]
        dip = np.load(f"{T}K/2_TDDFT/TRAJ/traj-{traj}/PLOTS_DATA/DIPOLE_RPA.dat.npy")[0,:,:].dot(pol)
        DIPOLE[Ti,traj] = dip[iSTATE]

E_AVE = np.average( ENERGY, axis=(1) )
print(E_AVE*27.2114)
# exit()
WC    = E_AVE
E_POL = np.zeros( (len(TEMP_LIST), NAVERAGES, NPOL) )
IPR   = np.zeros( (len(TEMP_LIST), NAVERAGES, NPOL) )
PHOT  = np.zeros( (len(TEMP_LIST), NAVERAGES, NPOL) )
for Ti,T in enumerate( TEMP_LIST ):
    for avei in range( NAVERAGES ):
        H = np.zeros( (NPOL,NPOL) )
        for A in range( NMOL ):
            H[A,A]  = ENERGY[Ti,avei+A]
            H[A,-1] = np.sqrt(WC[Ti]/2) * lambda_single * DIPOLE[Ti,avei+A]
            H[-1,A] = H[A,-1]
        H[-1,-1] = WC[Ti]
        
        e,u              = np.linalg.eigh( H )
        E_POL[Ti,avei,:] = e
        NORM             = np.sqrt( np.sum( u[:-1,:]**2, axis=0 ) )
        PROB             = u[:-1,:]**2 / NORM[None,:]
        IPR[Ti,avei,:]   = 1 / np.sum( PROB**2, axis=0 )
        PHOT[Ti,avei,:]  = u[-1,:]**2

for Ti,T in enumerate( TEMP_LIST ):
    plt.scatter( E_POL[Ti,:,:].flatten()*27.2114, IPR[Ti,:,:].flatten(), s=1, alpha=0.1, label="T = %d K" % (T) )
plt.xlim(6.4,7.4)
plt.legend()
plt.xlabel("Energy (eV)",fontsize=15)
plt.ylabel("IPR",fontsize=15)
plt.savefig("IPR.jpg", dpi=300)
plt.clf()
plt.close()

for Ti,T in enumerate( TEMP_LIST ):
    plt.scatter( E_POL[Ti,:,-1].flatten()*27.2114, IPR[Ti,:,-1].flatten(), s=1, alpha=0.1, label="T = %d K" % (T) )
plt.xlim(6.4,7.4)
plt.legend()
plt.xlabel("Energy (eV)",fontsize=15)
plt.ylabel("IPR",fontsize=15)
plt.savefig("IPR_UP.jpg", dpi=300)
plt.clf()
plt.close()

for Ti,T in enumerate( TEMP_LIST ):
    plt.scatter( E_POL[Ti,:,:].flatten()*27.2114, PHOT[Ti,:,:].flatten(), s=1, alpha=0.1, label=f"T = {T} K" )
plt.legend()
plt.xlim(6.4,7.4)
plt.xlabel("Energy (eV)",fontsize=15)
plt.ylabel("Transmission",fontsize=15)
plt.savefig("PHOT.jpg", dpi=300)
plt.clf()
plt.close()

