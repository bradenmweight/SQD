# List of Parameters for SQD Functionality

## Required Keywords
| Keyword | Type | Default | Function |
| ------ | ------ | ------ | ------ |
| NSTATES | int | None | Defines the number of electronic states used for the calculation. For example, "1" paried with "ISTATE = 0" requests Born-Oppenheimer molecular dynamics in the ground state. (TODO: Add functionality for non-GS BOMD.) |
| FUNCTIONAL | str | None | Defines the exchange-correlation functional for electronic structure. Must be understood by Gaussian16, *e.g.*, "B3LYP" |
| BASIS | str | None | Defines the atom-centered basis set for electronic structure. Must be understood by Gaussian16, *e.g.*, "6-31G*". (TODO: Add functionality for mixed basis as C:3-1G/Fe:LANL2DZ) |
| NSTEPS | int | None | Defines the number of nuclear time-steps for the molecular dynamics. |
| dtI | float | None | Defined the nuclear time-step in molecular dynamics. [fs] |
| ESTEPS | int | None | Defines the number of electronic time-steps per nuclear time-step for the quantum dynamics. |
| ISTATE | int | None | Defines the initial electronic state for the quantum dynamics [0,1,2,...] |
| CHARGE | int | None | Defines the total charge of the system. |
| MULTIPLICITY | int | None | Defines the electronic multiplicity of the system. |
| MEMORY | int | None | Defines the electronic multiplicity of the system. [GB] |
| NAMD_METHOD | str | None | Defines the semi-classical or mixed quantum-classical (MQC)quantum dynamics method for the propagation of the electronic and nuclear degrees of freedom. ["EH", "spinLSC", "GFSH"] See documentation for the description of each of thesde methods. |
| MD_ENSEMBLE | str | None | Defines the statistical mechanical ensemble for the molecular dyanmics. ["NVE", "NVT"] |

#
#
#
## Optional Keywords
| Keyword | Type | Default | Function |
| ------ | ------ | ------ | ------ |
| CPA | bool | False | Determines whether to do the classical path approximation (CPA) for the nuclear forces. |
