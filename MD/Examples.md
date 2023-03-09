# Examples of Various SQD Functionality

[DOCS]:   <https://bradenmweight.github.io/SQD/read.html?filename=Documentation.md>
[PARAMS]: <https://bradenmweight.github.io/SQD/read.html?filename=Parameters.md>

## Example Input Scripts

### BOMD in $S_0$
```
# Ground State (S0) BOMD with Langevin Bath (NVT)
NSTATES = 1              # Number of states in the calculation [int]
ISTATE  = 0              # Initial electronic state: 0, 1, 2, ... [int]
NSTEPS  = 10000          # Total MD steps [int]
dtI = 1.0                # Nuclear time-step (fs) [float]
FUNCTIONAL = AM1         # DFT functional or Semi-empirical method [str]
BASIS = STO-5G           # Local basis set [str]
CHARGE = 0               # System's net charge [int]
MULTIPLICITY = 1         # System electronic multiplicity [int]
MEMORY = 5               # Memory for Gaussian (GB) [int]
ESTEPS = 0               # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE = NVT        # NVE - Controls thermodynamic ensemble [str]
VELOC = ZERO             # Initial values for the nuclear velocities [str]
NCPUS_G16 = 24           # Number of shared-memory cores for QM calculation [int]
NVT_TYPE = LANGEVIN      # Lengevin thermostat [str]
LANGEVIN_LAMBDA = 50.0   # Coupling strength for Langevin thermostat (meV) [float]
TEMP = 300               # Temperature of the thermal bath (K) [float]
```


