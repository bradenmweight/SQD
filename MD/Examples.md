# List of Parameters for SQD Functionality

[DOCS]:   <https://bradenmweight.github.io/SQD/read.html?filename=Documentation.md>
[PARAMS]: <https://bradenmweight.github.io/SQD/read.html?filename=Parameters.md>
## Example Input Scripts

### BOMD in $S_0$
```
# Ground State (S0) BOMD with Langevin Bath (NVT)
NSTATES = 1 ! Includes ground state. NSTATES = 1 for GS BOMD [int]
ISTATE  = 0 # Initial electronic state for NAMD: 0, 1, 2, ... [int]
NSTEPS  = 10000 # Total MD steps [int]
dtI = 1.0 # Nuclear time-step (fs) [float]
FUNCTIONAL = AM1 # wB97XD # SVWN # DFT functional [str]
BASIS = STO-5G # Local basis set [str]
CHARGE = 0 # System's net charge [int]
MULTIPLICITY = 1 # System electronic multiplicity [int]
MEMORY = 5 # Memory for Gaussian (GB) [int]
ESTEPS = 0 # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE = NVT # NVE - Controls thermodynamic ensemble [str]
VELOC = ZERO # {ZERO} {MB} {READ} -> read from "input_veloc.xyz", [str]
NCPUS_G16 = 24 # [int]
NVT_TYPE = LANGEVIN # LANGEVIN, RESCALE [str]
LANGEVIN_LAMBDA = 50.0 # meV [float]
#RESCALE_FREQ = 25 # Steps between rescaling for NVT_TYPE = RESCALE [int]
TEMP = 300 # K [float]
```


