# Examples of Various SQD Functionality

[DOCS]:   <https://bradenmweight.github.io/SQD/read.html?filename=Documentation.md>
[PARAMS]: <https://bradenmweight.github.io/SQD/read.html?filename=Parameters.md>

## Example Workflow for Excited State Non-adiabtic Dynamics
1. Simulate Ground State Born-Oppenheimer Dyanmics at Constant Temperature (NVT) [Example 1]
2. Sample Ground State Thermal Distribution (scripts/extract_NVT_inital_conditions.py)
3. Copy "TRAJ" directory into a new NAMD folder
4. Populate and run all trajectories with NAMD.in [Example 3] and submission script (scripts/run_jobs.sh)
5. Extract results (scripts/get_average_population*, scripts/track_angle_allTRAJ.py, scripts/track_dihedral_angle_allTRAJ.py)
6. Publish a paper

## Example Input Scripts

### Example 1: BOMD in the Ground State with Langevin Bath (NVT)
```
# Ground State (S0) BOMD with Langevin Bath (NVT)
NSTATES         = 1          # Number of states in the calculation [int]
BOMD            = True       # Whether to do Born-Oppenheimer dynamics on ISTATE
ISTATE          = 0          # Initial electronic state: 0, 1, 2, ... [int]
NSTEPS          = 10000      # Total MD steps [int]
dtI             = 1.0        # Nuclear time-step (fs) [float]
FUNCTIONAL      = AM1        # DFT functional or Semi-empirical method [str]
BASIS           = STO-5G     # Local basis set [str]
CHARGE          = 0          # System's net charge [int]
MULTIPLICITY    = 1          # System electronic multiplicity [int]
MEMORY          = 5          # Memory for Gaussian (GB) [int]
ESTEPS          = 0          # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE     = NVT        # NVE - Controls thermodynamic ensemble [str]
VELOC           = ZERO       # Initial values for the nuclear velocities [str]
NCPUS_G16       = 24         # Number of shared-memory cores for QM calculation [int]
NVT_TYPE        = LANGEVIN   # Lengevin thermostat [str]
LANGEVIN_LAMBDA = 50.0       # Coupling strength for Langevin thermostat (meV) [float]
TEMP            = 300        # Temperature of the thermal bath (K) [float]
```

#
#
#


### Example 2: BOMD in the First Excited State at Constant Energy (NVE)
```
# Excited State (S1) BOMD at constant energy (NVE)
NSTATES         = 1          # Number of states in the calculation [int]
BOMD            = True       # Whether to do Born-Oppenheimer dynamics on ISTATE
ISTATE          = 0          # Initial electronic state: 0, 1, 2, ... [int]
NSTEPS          = 10000      # Total MD steps [int]
dtI             = 1.0        # Nuclear time-step (fs) [float]
FUNCTIONAL      = AM1        # DFT functional or Semi-empirical method [str]
BASIS           = STO-5G     # Local basis set [str]
CHARGE          = 0          # System's net charge [int]
MULTIPLICITY    = 1          # System electronic multiplicity [int]
MEMORY          = 5          # Memory for Gaussian (GB) [int]
ESTEPS          = 0          # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE     = NVE        # NVE - Controls thermodynamic ensemble [str]
VELOC           = READ       # Initial values for the nuclear velocities [str]
                             # Read from "./velocity_input.xyz"
NCPUS_G16       = 24         # Number of shared-memory cores for QM calculation [int]
```

### Example 3: NAMD starting in excited state = 2 and Including States = [0,1,2,3,4] in the Dynamics
```
# Non-adiabtic Dynamics (NAMD) starting in (S2) and including S0-S4
NSTATES         = 5          # Number of states in the calculation [int]
BOMD            = False      # Whether to do Born-Oppenheimer dynamics on ISTATE
ISTATE          = 2          # Initial electronic state: 0, 1, 2, ... [int]
NSTEPS          = 1000       # Total MD steps [int]
dtI             = 0.25       # Nuclear time-step (fs) [float]
FUNCTIONAL      = AM1        # DFT functional or Semi-empirical method [str]
BASIS           = STO-5G     # Local basis set [str]
CHARGE          = 0          # System's net charge [int]
MULTIPLICITY    = 1          # System electronic multiplicity [int]
MEMORY          = 5          # Memory for Gaussian (GB) [int]
NAMD_METHOD     = EH         # Controls which semi-classical/MQC method to use [str]
ESTEPS          = 200        # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE     = NVE        # NVE - Controls thermodynamic ensemble [str]
VELOC           = READ       # Initial values for the nuclear velocities [str]
                             # Read from "./velocity_input.xyz"
NCPUS_G16       = 1          # Number of shared-memory cores for QM calculation [int]
NCPUS_NAMD      = 3          # NCPUS, Dictates parallization across excited state forces. [int]
EL_PROP         = VV         # RK (4th-order explicit Runge-Kutta), VV (second-order symplectic) [str]
```

### Example 4: NAMD (in Classical Path Approximation) starting in excited state = 2 and Including States = [0,1,2,3,4] in the Dynamics
```
# Non-adiabtic Dynamics (NAMD) starting in (S2) and including S0-S4
NSTATES         = 5          # Number of states in the calculation [int]
BOMD            = False      # Whether to do Born-Oppenheimer dynamics on ISTATE
CPA             = True       # Whether to do Classical Path Approximation
ISTATE          = 2          # Initial electronic state: 0, 1, 2, ... [int]
NSTEPS          = 1000       # Total MD steps [int]
dtI             = 0.25       # Nuclear time-step (fs) [float]
FUNCTIONAL      = AM1        # DFT functional or Semi-empirical method [str]
BASIS           = STO-5G     # Local basis set [str]
CHARGE          = 0          # System's net charge [int]
MULTIPLICITY    = 1          # System electronic multiplicity [int]
MEMORY          = 5          # Memory for Gaussian (GB) [int]
NAMD_METHOD     = EH         # Controls which semi-classical/MQC method to use [str]
ESTEPS          = 200        # Number of electronic time steps per nuclear step [int]
MD_ENSEMBLE     = NVE        # NVE - Controls thermodynamic ensemble [str]
VELOC           = READ       # Initial values for the nuclear velocities [str]
                             # Read from "./velocity_input.xyz"
NCPUS_G16       = 1          # Number of shared-memory cores for QM calculation [int]
NCPUS_NAMD      = 3          # NCPUS, Dictates parallization across excited state forces. [int]
EL_PROP         = VV         # RK (4th-order explicit Runge-Kutta), VV (second-order symplectic) [str]
```