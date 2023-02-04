# here is a highlight
This is a test. Blah Blah


# Usage  
### Step 1
Create a folder and git clone this repository.
```
git clone https://github.com/arkajitmandal/SplitOperator1D 
```

### Step 2
Code up the model system in a python file inside the "Model" folder and name it  'modelName.py'.  

The 'modelName.py' should look like:
```py
import numpy as np
from tools import Diag 
from numpy import kron as ꕕ


class parameters():
    Rmin  = -13
    Rmax  = 13
    nR    = 512
    dt    = 1.0
    steps = 1200
    aniskip = 250
    M = 2000
#--------------------------------------------------------
#-------------------   Model        ---------------------
#--------------------------------------------------------

def Hel(R):
    """Hel for Tully 1

    Args:
        R (float): nuclear position

    Returns:
        N x N matrix: Matrix elements of electronic part of the Hamiltonian
    """


    V = np.zeros((2,2))
    A = 0.01
    B = 1.6
    C = 0.005
    D = 1.0
    V[0,0] = A*np.sign(R)*(1-np.exp(-B*np.abs(R)))
    V[1,1] = -V[0,0]
    V[0,1] = C*np.exp(-D*R**2)
    V[1,0] = V[0,1]

    return V

#--------------------------------------------------------
#-------------------   initial Ψ    ---------------------
#--------------------------------------------------------
def psi(R):
  """Initial wavefunction

  Args:
      R (numpy array): R is a numerical grid over which the nuclear part of
      the wavefuntion is evaluated. 

  Returns:
      Ψ: wavefunction in the nuclear ⊗ electronic wavefunction. 
      I have used a initial state: Ψ(R) = χ(R) ⊗ |i><i| = χ(R) ⊗ φ
      can be easily modified to have a entangled state--> Ψ(R) = ∑ χi(R) ⊗ |i><i|
  """
  
  # Nuclear Part
  α = 1.0 
  R0 =  -9.0
  P0 = 30
  χ = np.exp(- 0.5 * α * (R - R0)**2.0 ) * np.exp( 1j * P0 * (R - R0))
  χ =  χ/(np.sum(χ*χ.conjugate())**0.5)
  
  # Electronic Part
  Φ = np.array([1, 0])
  return  ꕕ(Φ, χ)
```
You can find several examples of model files inside the "Model" folder. I will explain each parts of this file in more detain in a section below.
