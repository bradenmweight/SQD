import numpy as np

### NAMD METHODS ###
import Eh
import spinLSC
import GFSH

def initialize_mapping(DYN_PROPERTIES):
    """
    Wrapper for initializing mapping variables
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.initialize_mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.initialize_mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "GFSH" ):
        return GFSH.initialize_mapping(DYN_PROPERTIES)

def propagage_Mapping(DYN_PROPERTIES):
    """
    Wrapper for electronic propagation
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.propagage_Mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.propagage_Mapping(DYN_PROPERTIES)
    elif ( DYN_PROPERTIES["NAMD_METHOD"] == "GFSH" ):
        return GFSH.propagage_Mapping(DYN_PROPERTIES)

def rotate_Mapping(DYN_PROPERTIES):
    """
    Wrapper for the transformation of mapping variables
    """
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "EH" ):
        return Eh.rotate_Mapping(DYN_PROPERTIES)
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "SPINLSC" ):
        return spinLSC.rotate_Mapping(DYN_PROPERTIES)
    if ( DYN_PROPERTIES["NAMD_METHOD"] == "GFSH" ):
        return GFSH.rotate_Mapping(DYN_PROPERTIES)