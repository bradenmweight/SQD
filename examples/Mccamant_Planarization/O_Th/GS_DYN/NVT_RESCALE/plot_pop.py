import numpy as np
from matplotlib import pyplot as plt

tmp  = np.loadtxt("MD_OUTPUT/mapping_re.dat") + 1.0j * np.loadtxt("MD_OUTPUT/mapping_im.dat")
time = np.real(tmp[:,0])
z    = tmp[:,1:]
NStates = len(z[0,:])
rho  = np.einsum( "tj,tk->tjk", np.conjugate(z), z )


#for state in range(NStates):
    #plt.plot( time, np.real(rho[:,state,state]), label=f"S{state}" )
plt.plot( time, np.real(rho[:,0,0] + rho[:,1,1]), label=f"Total" )
plt.legend()
plt.savefig("MD_OUTPUT/POP.jpg", dpi=600)

