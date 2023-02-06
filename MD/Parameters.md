# List of Parameters for SQD Functionality

[DOCS]: <https://bradenmweight.github.io/SQD/read.html?filename=Documentation.md>

## Required Keywords
| Keyword | Type | Default | Function |
| ------ | ------ | ------ | ------ |
| NSTATES | int | None | Defines the number of electronic states used for the calculation. For example, "NSTATES = 1" paried with "ISTATE = 0" requests Born-Oppenheimer molecular dynamics (BOMD) in the ground state. (TODO: Add functionality for non-GS BOMD.) |
| ISTATE | int | None | Defines the initial electronic state for the quantum dynamics [0,1,2,...] |
| NSTEPS | int | None | Defines the number of nuclear time-steps for the molecular dynamics. |
| dtI | float | None | Defined the nuclear time-step in molecular dynamics. [fs] |
| ESTEPS | int | None | Defines the number of electronic time-steps per nuclear time-step for the quantum dynamics. |
| FUNCTIONAL | str | None | Defines the exchange-correlation functional for electronic structure. Must be understood by Gaussian16, *e.g.*, "B3LYP" |
| BASIS | str | None | Defines the atom-centered basis set for electronic structure. Must be understood by Gaussian16, *e.g.*, "6-31G*". (TODO: Add functionality for mixed basis as C:3-1G/Fe:LANL2DZ) |
| CHARGE | int | None | Defines the total charge of the system. |
| MULTIPLICITY | int | None | Defines the electronic multiplicity of the system. |
| MEMORY | int | None | Defines the amount of memory given to Gaussian16 for electronic structure calculations. [GB] |
| NAMD_METHOD | str | None | Defines the semi-classical or mixed quantum-classical (MQC) quantum dynamics method for the propagation of the electronic and nuclear degrees of freedom. ["EH", "spinLSC", "GFSH"] See [documentation][DOCS] for the description of each of thesde methods. |
| MD_ENSEMBLE | str | None | Defines the statistical mechanical ensemble for the molecular dyanmics. ["NVE", "NVT"] |



#
#
#
## Optional Keywords
| Keyword | Type | Default | Function |
| ------ | ------ | ------ | ------ |
| CPA | bool | False | Determines whether to do the classical path approximation (CPA) for the nuclear forces. |
| PARALLEL_FORCES | bool | False | Determines whether to parallelize (across excited state) for the calculation of nuclear forces. For example, if there are 5 electronic excited state in the calculation, with "PARALLEL_FORCES = True" and "NCPUS_NAMD  = 5", then SQD will simultaneously run all excited state calculations in parallel, else they will be run in serial.
| NCPUS_NAMD | int | None | Determines the number of CPUs available for excited state parallelization. See "PARALLEL_FORCES". |
| NCPUS_G16 | int | 1 | Determines the number of CPUs available for Gaussian16 electronic structure calculations. |
| RUN_ELEC_STRUC | str | USE_CURRENT_NODE | Determines how to submit Gaussian16 electronic structure calculations. ["USE_CURRENT_NODE","SUBMIT_SBATCH"] |
| SBATCH_G16 | str | None | Absolute file path to Gaussian16 submission script. "RUN_ELEC_STRUC" must be set to "SUBMIT_SBATCH". |
| EL_PROP | str | VV | Determines the time-integrator for the electronic degrees of freedom. ["VV","RK"] See [documentation][DOCS] for the description of each of these methods. |
| NVT_TYPE | str | LANGEVIN | Determines the type of thermal bath in canonical ensemble (*i.e.*, NVT) molecular dynamics. ["LANGEVIN","RESCALE"] See [documentation][DOCS] for the description of each of these methods. |
| LANGEVIN_LAMBDA | float | None | Determines the coupling/friction parameter in classical Langevin molecular dynamics. "NVT_TYPE" must be "LANGEVIN". See [documentation][DOCS] for the description of each of these methods. |
| TEMP | float | None | Determines the target temperature in canonical ensemble (*i.e.*, NVT) molecular dynamics. See [documentation][DOCS] for the description of each of these methods. |
| RESCALE_FREQ | int | None | Determines the frequency in which to rescale the nuclear velocities to achieve the target temperature in canonical ensemble (*i.e.*, NVT) molecular dynamics. "NVT_TYPE" must be "RESCALE". See [documentation][DOCS] for the description of each of these methods. |
| VELOC | str | MB | Determines how the classical nuclear velocities are initialized. ["ZERO","MB","READ"] See [documentation][DOCS] for the description of each of these methods. |
| DATA_SAVE_FREQ | int | 1 | Determines how often the main data are saved to text files. |


