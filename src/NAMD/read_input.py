import numpy as np
import subprocess as sp
import os

def read():

    DYN_PROPERTIES = {} # Everything will come from here

    # Read input file
    input_lines = open('NAMD.in','r').readlines()
    for count, line in enumerate(input_lines):
        ### Clean line and Check for comments ###
        t = line.split()
        if ( len(t) == 0 or line.split()[0] in ["#","!"] ): continue # Check for comment line
        t = [ j.strip() for j in line.split("=") ]
        tnew = []
        for tj in t:
            if ( "#" in tj.split() ):
                tnew.append( tj.split()[:tj.split().index("#")][0] )
                break
            if ( "!" in tj.split() ):
                tnew.append( tj.split()[:tj.split().index("!")][0] )
                break
            else:
                tnew.append(tj)
        t = tnew

        if ( len(t) == 2 ):
            
            # Look for NStates
            if ( t[0].upper() == "nstates".upper() ):
                try:
                    DYN_PROPERTIES["NStates"] = int( t[1] )
                except ValueError:
                    print(f"\t'NSteps' must be an integer: '{t[1]}'")
                    exit()

            # Look for NSteps
            if ( t[0].upper() == "nsteps".upper() ):
                try:
                    DYN_PROPERTIES["NSteps"] = int( t[1] )
                except ValueError:
                    print(f"\t'NSteps' must be an integer: '{t[1]}'")
                    exit()

            # Look for ISTATE
            if ( t[0].upper() == "istate".upper() ):
                try:
                    DYN_PROPERTIES["ISTATE"] = int( t[1] )
                except ValueError:
                    print(f"\t'ISTATE' must be an integer: '{t[1]}'")
                    exit()

            # Look for FUNCTIONAL
            if ( t[0].upper() == "FUNCTIONAL".upper() ):
                DYN_PROPERTIES["FUNCTIONAL"] = t[1].upper()
                # The accuracy of input is up to the user...scary!

            # Look for BASIS
            if ( t[0].upper() == "BASIS".upper() ):
                DYN_PROPERTIES["BASIS_SET"] = t[1].upper()
                # The accuracy of input is up to the user...scary!
            # TODO -- Add functionality for mixed basis sets

            # Look for ISTATE
            if ( t[0].upper() == "dtI".upper() ):
                try:
                    DYN_PROPERTIES["dtI"] = float( t[1] ) * 41.341 # 41.341 a.u. / fs
                except ValueError:
                    print(f"\t'dtI' must be a float: '{t[1]}'")
                    exit()
            
            # Look for ESTEPS
            if ( t[0].upper() == "ESTEPS".upper() ):
                try:
                    DYN_PROPERTIES["ESTEPS"] = int( t[1] )
                except ValueError:
                    print(f"\t'ESTEPS' must be an integer: '{t[1]}'")
                    exit()

            # Look for NCPUS
            if ( t[0].upper() == "NCPUS_NAMD".upper() ):
                try:
                    DYN_PROPERTIES["NCPUS_NAMD"] = int( t[1] )
                except ValueError:
                    print(f"\t'NCPUS' must be an integer: '{t[1]}'")
                    exit()

            # Look for NCPUS
            if ( t[0].upper() == "NCPUS_G16".upper() ):
                try:
                    DYN_PROPERTIES["NCPUS_G16"] = int( t[1] )
                except ValueError:
                    print(f"\t'NCPUS' must be an integer: '{t[1]}'")
                    exit()

            # Look for MEMORY
            if ( t[0].upper() == "MEMORY".upper() ):
                try:
                    DYN_PROPERTIES["MEMORY"] = int( t[1] )
                except ValueError:
                    print(f"\t'MEMORY' must be an integer (units of GB): '{t[1]}'")
                    exit()

            # Look for CHARGE
            if ( t[0].upper() == "CHARGE".upper() ):
                try:
                    DYN_PROPERTIES["CHARGE"] = int( t[1] )
                except ValueError:
                    print(f"\t'CHARGE' must be an integer: '{t[1]}'")
                    exit()

            # Look for CHARGE
            if ( t[0].upper() == "MULTIPLICITY".upper() ):
                try:
                    DYN_PROPERTIES["MULTIPLICITY"] = int( t[1] )
                except ValueError:
                    print(f"\t'MULTIPLICITY' must be an integer: '{t[1]}'")
                    exit()

            # Look for VELOC
            if ( t[0].upper() == "VELOC".upper() ):
                DYN_PROPERTIES["VELOC"] = t[1].upper()
                # Later, we will check this input.

            # Look for RUN_ELEC_STRUC
            if ( t[0].upper() == "RUN_ELEC_STRUC".upper() ):
                DYN_PROPERTIES["RUN_ELEC_STRUC"] = t[1].upper()
                # The accuracy of input is up to the user...scary!

            # Look for SBATCH_G16
            if ( t[0].upper() == "SBATCH_G16".upper() ):
                DYN_PROPERTIES["SBATCH_G16"] = t[1] # DO NOT CHANGE CASE
                # The accuracy of input is up to the user...scary!

            # Look for PARALLEL_FORCES
            if ( t[0].upper() == "PARALLEL_FORCES".upper() ):
                try:
                    DYN_PROPERTIES["PARALLEL_FORCES"] = bool( t[1] )
                except ValueError:
                    print("Input for 'PARALLEL_FORCES' must be a boolean. (True or False)")
                    exit()

            # Look for FORCE_MAP_NORM -- FOR DEBUGGING ONLY
            if ( t[0].upper() == "FORCE_MAP_NORM".upper() ):
                try:
                    DYN_PROPERTIES["FORCE_MAP_NORM"] = bool( t[1] )
                except ValueError:
                    print("Input for 'FORCE_MAP_NORM' must be a boolean. (True or False)")
                    exit()

            # Look for EL_PROP
            if ( t[0].upper() == "EL_PROP".upper() ):
                DYN_PROPERTIES["EL_PROP"] = t[1].upper()
                if ( DYN_PROPERTIES["EL_PROP"] not in ["VV","RK"] ):
                    print("Input for 'EL_PROP' must be either 'VV' (Velocity-Verlet) or 'RK' (Runge-Kutta).")
                    exit()

            # Look for NAMD_METHOD
            if ( t[0].upper() == "NAMD_METHOD".upper() ):
                DYN_PROPERTIES["NAMD_METHOD"] = t[1].upper()
                if ( DYN_PROPERTIES["NAMD_METHOD"] not in ["EH","SPINLSC","GFSH"] ):
                    print("Input for 'NAMD_METHOD' must be either 'EH', 'SPINLSC', or 'GFSH'.")
                    exit()

            # Look for MD_ENSEMBLE
            if ( t[0].upper() == "MD_ENSEMBLE".upper() ):
                DYN_PROPERTIES["MD_ENSEMBLE"] = t[1].upper()
                if ( DYN_PROPERTIES["MD_ENSEMBLE"] not in ["NVT","NVE"] ):
                    print("Input for 'MD_ENSEMBLE' must be either 'NVT' or 'NVE'.")
                    exit()

            # Look for NVT_TYPE
            if ( t[0].upper() == "NVT_TYPE".upper() ):
                DYN_PROPERTIES["NVT_TYPE"] = t[1].upper()
                if ( DYN_PROPERTIES["NVT_TYPE"] not in ["LANGEVIN", "RESCALE"] ):
                    print("Input for 'NVT_TYPE' must be 'LANGEVIN' or 'RESCALE'.")
                    exit()

            # Look for LANGEVIN_LAMBDA
            if ( t[0].upper() == "LANGEVIN_LAMBDA".upper() ):
                try:
                    DYN_PROPERTIES["LANGEVIN_LAMBDA"] = float( t[1] )
                except ValueError:
                    print(f"\t'LANGEVIN_LAMBDA' must be a float: '{t[1]}'")
                    exit()
                assert( DYN_PROPERTIES["LANGEVIN_LAMBDA"] >= 0 ), f"'LANGEVIN_LAMBDA' must be greater than or equal to 0.0: {DYN_PROPERTIES['LANGEVIN_LAMBDA']}"

            # Look for RESCALE_FREQ
            if ( t[0].upper() == "RESCALE_FREQ".upper() ):
                try:
                    DYN_PROPERTIES["RESCALE_FREQ"] = int( t[1] )
                except ValueError:
                    print(f"\t'RESCALE_FREQ' must be an integer: '{t[1]}'")
                    exit()
                assert( DYN_PROPERTIES["RESCALE_FREQ"] >= 0 ), f"'RESCALE_FREQ' must be greater than or equal to 0: {DYN_PROPERTIES['RESCALE_FREQ']}"


            # Look for TEMP
            if ( t[0].upper() == "TEMP".upper() ):
                try:
                    DYN_PROPERTIES["TEMP"] = float( t[1] )
                except ValueError:
                    print(f"\t'TEMP' must be a float: '{t[1]}'")
                    exit()
                assert( DYN_PROPERTIES["TEMP"] >= 0 ), f"'TEMP' must be greater than or equal to 0.0: {DYN_PROPERTIES['TEMP']}"


            # Look for DATA_SAVE_FREQ
            if ( t[0].upper() == "DATA_SAVE_FREQ".upper() ):
                try:
                    DYN_PROPERTIES["DATA_SAVE_FREQ"] = int( t[1] )
                except ValueError:
                    print(f"\t'DATA_SAVE_FREQ' must be an integer: '{t[1]}'")
                    exit()


            # Look for EL_INTERPOLATION
            if ( t[0].upper() == "EL_INTERPOLATION".upper() ):
                try:
                    DYN_PROPERTIES["EL_INTERPOLATION"] = bool( t[1] )
                except ValueError:
                    print(f"\t'EL_INTERPOLATION' must be a bool: '{t[1]}'")
                    exit()

            # Look for TDDFT_CONVERG
            if ( t[0].upper() == "TDDFT_CONVERG".upper() ):
                try:
                    DYN_PROPERTIES["TDDFT_CONVERG"] = int( t[1] )
                except ValueError:
                    print(f"\t'TDDFT_CONVERG' must be an int: '{t[1]}'")
                    exit()


            # Look for CHECK_TRIVIAL_CROSSING
            if ( t[0].upper() == "CHECK_TRIVIAL_CROSSING".upper() ):
                try:
                    DYN_PROPERTIES["CHECK_TRIVIAL_CROSSING"] = bool( t[1] )
                except ValueError:
                    print(f"\t'CHECK_TRIVIAL_CROSSING' must be a bool: '{t[1]}'")
                    exit()

            # Look for PRINT_NACR
            if ( t[0].upper() == "PRINT_NACR".upper() ):
                try:
                    DYN_PROPERTIES["PRINT_NACR"] = bool( t[1] )
                except ValueError:
                    print(f"\t'PRINT_NACR' must be a bool: '{t[1]}'")
                    exit()

            # Look for BOMD
            if ( t[0].upper() == "BOMD".upper() ):
                if ( t[1].upper() == "TRUE" ):
                    DYN_PROPERTIES["BOMD"] = True
                elif ( t[1].upper() == "FALSE" ):
                    DYN_PROPERTIES["BOMD"] = False
                else:
                    print("Input for 'BOMD' must be a boolean. (True or False)")
                    exit()

            # Look for SYMMETRIZE_OVERLAP
            if ( t[0].upper() == "SYMMETRIZE_OVERLAP".upper() ):
                if ( t[1].upper() == "TRUE" ):
                    DYN_PROPERTIES["SYMMETRIZE_OVERLAP"] = True
                elif ( t[1].upper() == "FALSE" ):
                    DYN_PROPERTIES["SYMMETRIZE_OVERLAP"] = False
                else:
                    print("Input for 'SYMMETRIZE_OVERLAP' must be a boolean. (True or False)")
                    exit()

            # Look for S0S1_CI_FLAG
            if ( t[0].upper() == "S0S1_CI_FLAG".upper() ):
                if ( t[1].upper() == "TRUE" ):
                    DYN_PROPERTIES["S0S1_CI_FLAG"] = True
                elif ( t[1].upper() == "FALSE" ):
                    DYN_PROPERTIES["S0S1_CI_FLAG"] = False
                else:
                    print("Input for 'S0S1_CI_FLAG' must be a boolean. (True or False)")
                    exit()


            # Look for CPA
            if ( t[0].upper() == "CPA".upper() ):
                try:
                    if ( t[1].upper() == "TRUE" ):
                        DYN_PROPERTIES["CPA"] = True
                    elif ( t[1].upper() == "FALSE" ):
                        DYN_PROPERTIES["CPA"] = False
                    else:
                        print("Input for 'CPA' must be a boolean. (True or False)")
                        exit()
                except:
                    print("Something wrong with CPA block in 'read_input.py'")
                    exit()
        else:
            print( f"Error: Input is wrong at line {count+1}: {line}" )
            print( f"\tToo many '='" )
            exit()

    try:
        print("MANDATORY INPUT VARIABLES:")
        print( "  NStates ="); print("\t\t", DYN_PROPERTIES["NStates"] )
        print( "  NSteps ="); print("\t\t", DYN_PROPERTIES["NSteps"] )
        print( "  dtI ="); print("\t\t", DYN_PROPERTIES["dtI"]/41.341, "(fs)" )
        print( "  ESTEPS ="); print("\t\t", DYN_PROPERTIES["ESTEPS"] )
        print( "  ISTATE ="); print("\t\t", DYN_PROPERTIES["ISTATE"] )
        print( "  FUNCTIONAL ="); print("\t\t", DYN_PROPERTIES["FUNCTIONAL"] )
        if( DYN_PROPERTIES["FUNCTIONAL"] not in ["DFTB", "DFTBA"] ): print( "  BASIS_SET ="); print("\t\t", DYN_PROPERTIES["BASIS_SET"] )
        print( "  CHARGE ="); print("\t\t", DYN_PROPERTIES["CHARGE"] )
        print( "  MULTIPLICITY ="); print("\t\t", DYN_PROPERTIES["MULTIPLICITY"] )
        print( "  MEMORY ="); print("\t\t", DYN_PROPERTIES["MEMORY"], "(GB)" )
        print( "  NAMD_METHOD ="); print("\t\t", DYN_PROPERTIES["NAMD_METHOD"] )
        print( "  EL_PROP ="); print("\t\t", DYN_PROPERTIES["EL_PROP"] )
        print( "  MD_ENSEMBLE ="); print("\t\t", DYN_PROPERTIES["MD_ENSEMBLE"] )
        print( "  VELOC ="); print("\t\t", DYN_PROPERTIES["VELOC"] )
    except KeyError:
        print("Input file is missing a mandatory entry (see above). Check it.")
        exit()

    #assert( DYN_PROPERTIES["ISTATE"] <= DYN_PROPERTIES["NStates"]-1 ), "ISTATE must be less than the total number of states."

    # Try to read the SQD_HOME path
    SQD_PATH = sp.check_output("echo $SQD_HOME",shell=True).decode().strip("\n")
    if ( SQD_PATH == "" ):
        print("SQD home path is not set.\n\texport SQD_HOME=/absolute/path/to/SQD/")
        exit()
    else:
        DYN_PROPERTIES["SQD_HOME_PATH"] = SQD_PATH

    # Try to read the SQD_SCRATCH path
    SQD_SCRATCH = sp.check_output("echo $SQD_SCRATCH",shell=True).decode().strip("\n")
    if ( SQD_SCRATCH == "" ):
        print("SQD scratch path is not set.\n\texport SQD_SCRATCH=/absolute/path/to/SQD_SCRATCH/")
        exit()
    else:
        DYN_PROPERTIES["SQD_SCRATCH_PATH"] = SQD_SCRATCH
    
    # Save running directory
    DYN_PROPERTIES["SQD_RUNNING_DIR"] = os.getcwd()
    
    
    
    return DYN_PROPERTIES



def read_geom():
    XYZ_File = open("geometry_input.xyz","r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_labels = []
    Atom_coords_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_labels.append( t[0] )
        Atom_coords_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ]) / 0.529 # Ang -> a.u.

    assert( int(XYZ_File[0]) == len(Atom_labels) ), "Number of atoms incorrect in input XYZ file."
    assert( int(XYZ_File[0]) == len(Atom_coords_new) ), "Number of atoms incorrect in input XYZ file."

    return Atom_labels, Atom_coords_new

def read_veloc():
    """
    TODO Add checks for XYZ user input
    """
    XYZ_File = open("velocity_input.xyz","r").readlines()
    NAtoms = int(XYZ_File[0])
    Atom_velocs_new = np.zeros(( NAtoms, 3 ))
    for count, line in enumerate(XYZ_File[2:]):
        t = line.split()
        Atom_velocs_new[count,:] = np.array([ float(t[1]), float(t[2]), float(t[3]) ]) / 0.529 / 41.341 # Ang -> a.u.

    return Atom_velocs_new

def set_masses(Atom_labels):
    mass_amu_to_au = 1837/1.007 # au / amu
    masses_amu = \
{"H":   1.00797,
"He":	4.00260,
"Li":	6.941,
"Be":	9.01218,
"B":    10.81,
"C":    12.011,
"N":    14.0067,
"O":    15.9994,
"F":    18.998403,
"Ne":	20.179,
"Na":	22.98977,
"Mg":	24.305,
"Al":	26.98154,
"Si":	28.0855,
"P":    30.97376,
"S":    32.06,
"Cl":	35.453,
"K":    39.0983,
"Ar":	39.948,
"Ca":	40.08,
"Sc":	44.9559,
"Ti":	47.90,
"V":    50.9415,
"Cr":	51.996,
"Mn":	54.9380,
"Fe":	55.847,
"Ni":	58.70,
"Co":	58.9332,
"Cu":	63.546,
"Zn":	65.38,
"Ga":	69.72,
"Ge":	72.59,
"As":	74.9216,
"Se":	78.96,
"Br":	79.904,
"Kr":	83.80,
"Rb":	85.4678,
"Sr":	87.62,
"Y":    88.9059,
"Zr":	91.22,
"Nb":	92.9064,
"Mo":	95.94,
"Ru":	101.07,
"Rh":	102.9055,
"Pd":	106.4,
"Ag":	107.868,
"Cd":	112.41,
"In":	114.82,
"Sn":	18.69,
"Sb":	121.75,
"I":    126.9045,
"Te":	127.60,
"Xe":	131.30,
"Cs":	132.9054,
"Ba":	137.33,
"La":	138.9055,
"Ce":	140.12,
"Pr":	140.9077,
"Nd":	144.24,
"Sm":	150.4,
"Eu":	151.96,
"Gd":	157.25,
"Tb":	158.9254,
"Dy":	162.50,
"Ho":	164.9304,
"Er":	167.26,
"Tm":	168.9342,
"Yb":	173.04,
"Lu":	174.967,
"Hf":	178.49,
"Ta":	180.9479,
"W":    183.85,
"Re":	186.207,
"Os":	190.2,
"Ir":	192.22,
"Pt":	195.09,
"Au":	196.9665,
"Hg":	200.59,
"Tl":	204.37,
"Pb":	207.2,
"Bi":	208.9804,
"Ra":	226.0254}


    masses = []
    for at in Atom_labels:
        masses.append( masses_amu[at] )
    return np.array(masses) * mass_amu_to_au

def get_initial_velocs(DYN_PROPERTIES):
    
    Atom_labels = DYN_PROPERTIES["Atom_labels"]
    masses      = DYN_PROPERTIES["MASSES"]

    # TODO Get Wigner distribution for initial velocities

    if ( DYN_PROPERTIES["VELOC"] == "MB" ):
        import random
        velocs = np.zeros(( len(Atom_labels), 3 ))
        T = 300 # K
        kT  = T * (0.025/300) / 27.2114 # K -> eV -> au
        V0  = np.sqrt( kT / masses )
        SIG = kT / masses
        for at,atom in enumerate(Atom_labels):
            for d in range(3):
                velocs[at,d] = random.gauss( V0[at], SIG[at] )
    
    elif (DYN_PROPERTIES["VELOC"] == "ZERO"):
        velocs = np.zeros(( len(Atom_labels), 3 ))
    
    elif ( DYN_PROPERTIES["VELOC"] == "READ" ): # This will be usual way to perform with Wigner for now.
        velocs = read_veloc() # Reads "velocity_input.xyz"

    else:
        assert(False), "Initial velocities not specified properly.\t Must be 'ZERO', 'MB', or 'RED'."

    return velocs

def initialize_MD_variables(DYN_PROPERTIES):
    
    DYN_PROPERTIES["MD_STEP"] = 0    
    DYN_PROPERTIES["Atom_labels"], DYN_PROPERTIES["Atom_coords_new"] = read_geom()
    DYN_PROPERTIES["NAtoms"] = len( DYN_PROPERTIES["Atom_labels"] )
    DYN_PROPERTIES["MASSES"] = set_masses(DYN_PROPERTIES["Atom_labels"])
    DYN_PROPERTIES["dtE"]    = DYN_PROPERTIES["dtI"] / DYN_PROPERTIES["ESTEPS"]
    DYN_PROPERTIES["Atom_velocs_new"] = get_initial_velocs(DYN_PROPERTIES)

    

    try:
        tmp = DYN_PROPERTIES["RUN_ELEC_STRUC"]
    except KeyError:
        DYN_PROPERTIES["RUN_ELEC_STRUC"] = "use_current_node".upper()

    try:
        tmp = DYN_PROPERTIES["EL_PROP"]
    except KeyError:
        DYN_PROPERTIES["EL_PROP"] = "VV"

    try:
        tmp = DYN_PROPERTIES["NCPUS_NAMD"]
    except KeyError:
        DYN_PROPERTIES["NCPUS_NAMD"] = 1

    try:
        tmp = DYN_PROPERTIES["NCPUS_G16"]
    except KeyError:
        DYN_PROPERTIES["NCPUS_G16"] = 1

    try:
        tmp = DYN_PROPERTIES["REMOVE_COM_MOTION"]
    except KeyError:
        DYN_PROPERTIES["REMOVE_COM_MOTION"] = True # Default is to remove COM motion

    try:
        tmp = DYN_PROPERTIES["REMOVE_ANGULAR_VELOCITY"]
    except KeyError:
        DYN_PROPERTIES["REMOVE_ANGULAR_VELOCITY"] = True # Default is to remove angular velocity

    try:
        tmp = DYN_PROPERTIES["PARALLEL_FORCES"]
    except KeyError:
        DYN_PROPERTIES["PARALLEL_FORCES"] = False # Set to Flase by default

    try:
        tmp = DYN_PROPERTIES["PRINT_NACR"]
    except KeyError:
        DYN_PROPERTIES["PRINT_NACR"] = False # Set to Flase by default -- Large file

    try:
        tmp = DYN_PROPERTIES["TDDFT_CONVERG"]
    except KeyError:
        DYN_PROPERTIES["TDDFT_CONVERG"] = 4 # Default is to do 10**-N, N = 4

    if ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVT" ):
        
        try:
            tmp = DYN_PROPERTIES["NVT_TYPE"]
        except KeyError:
            assert(False), f"\t'NVT_TYPE' needs to be defined if 'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']}"
        
        try:
            tmp = DYN_PROPERTIES["LANGEVIN_LAMBDA"]
        except KeyError:
            if( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
                assert(False), f"\t'LANGEVIN_LAMBDA' needs to be defined if 'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']}"
            elif( DYN_PROPERTIES["NVT_TYPE"] == "RESCALE" ):
                DYN_PROPERTIES["LANGEVIN_LAMBDA"] = None
            else:
                assert(False), f"\t'MD_ENSEMBLE', 'NVT_TYPE', and 'RESCALE_FREQ' need to be consistent."
        
        try:
            tmp = DYN_PROPERTIES["RESCALE_FREQ"]
        except KeyError:
            if( DYN_PROPERTIES["NVT_TYPE"] == "LANGEVIN" ):
                DYN_PROPERTIES["RESCALE_FREQ"] = None
            elif( DYN_PROPERTIES["NVT_TYPE"] == "RESCALE" ):
                assert(False), f"\t'RESCALE_FREQ' needs to be defined if 'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']}"
            else:
                assert(False), f"\t'MD_ENSEMBLE', 'NVT_TYPE', and 'RESCALE_FREQ' need to match."

        try:
            tmp = DYN_PROPERTIES["TEMP"]
        except KeyError:
            assert(False), f"\t'TEMP' needs to be defined if 'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']}"
    
    elif ( DYN_PROPERTIES["MD_ENSEMBLE"] == "NVE" ):
        
        try:
            tmp = DYN_PROPERTIES["NVT_TYPE"]
        except KeyError:
            DYN_PROPERTIES["NVT_TYPE"] = None
        assert(DYN_PROPERTIES["NVT_TYPE"] == None ), f"\n\t'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']} and 'NVT_TYPE' = {DYN_PROPERTIES['NVT_TYPE']} are not compatible.\n\t 'NVT_TYPE' should not appear as we are trying to do NVE dynamics !"
        
        try:
            tmp = DYN_PROPERTIES["LANGEVIN_LAMBDA"]
        except KeyError:
            DYN_PROPERTIES["LANGEVIN_LAMBDA"] = None
        assert(DYN_PROPERTIES["LANGEVIN_LAMBDA"] == None ), f"\n\t'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']} and 'LANGEVIN_LAMBDA' = {DYN_PROPERTIES['LANGEVIN_LAMBDA']} are not compatible.\n\t'LANGEVIN_LAMBDA' should not appear as we are trying to do NVE dynamics !"
        
        try:
            tmp = DYN_PROPERTIES["TEMP"]
        except KeyError:
            DYN_PROPERTIES["TEMP"] = None
        assert(DYN_PROPERTIES["TEMP"] == None ), f"\n\t'MD_ENSEMBLE' = {DYN_PROPERTIES['MD_ENSEMBLE']} and 'TEMP' = {DYN_PROPERTIES['TEMP']} are not compatible.\n\t'TEMP' should not appear as we are trying to do NVE dynamics !"
    
    else:
        assert(False), "'MD_ENSEMBLE' needs to be defined as either 'NVE' or 'NVT'."
    


    try:
        tmp = DYN_PROPERTIES["DATA_SAVE_FREQ"]
    except KeyError:
        DYN_PROPERTIES["DATA_SAVE_FREQ"] = 1 # Default is to save every step. Might make large output files for NVT

    try:
        tmp = DYN_PROPERTIES["CPA"]
    except KeyError:
        DYN_PROPERTIES["CPA"] = False # Default is not to do classical path approximation

    try:
        tmp = DYN_PROPERTIES["BOMD"]
    except KeyError:
        DYN_PROPERTIES["BOMD"] = False # Default is not to do BOMD

    try:
        tmp = DYN_PROPERTIES["SYMMETRIZE_OVERLAP"]
    except KeyError:
        DYN_PROPERTIES["SYMMETRIZE_OVERLAP"] = True # Default is not to do

    try:
        tmp = DYN_PROPERTIES["S0S1_CI_FLAG"]
    except KeyError:
        DYN_PROPERTIES["S0S1_CI_FLAG"] = False # Default is not to do


    if ( DYN_PROPERTIES["CPA"] == True and DYN_PROPERTIES["BOMD"] == True ):
        print("\nError: CPA = True and BOMD = True")
        print("Error: This does not make sense.\n")
        exit()


    ########################################################################
    ############## BELOW ARE PRIMARILY FOR DEBUGGING PURPOSES ##############
    ########################################################################

    try:
        tmp = DYN_PROPERTIES["FORCE_MAP_NORM"]
    except KeyError:
        DYN_PROPERTIES["FORCE_MAP_NORM"] = False # Default is not to do this. This is for de-bugging purposes.

    try:
        tmp = DYN_PROPERTIES["EL_INTERPOLATION"]
    except KeyError:
        DYN_PROPERTIES["EL_INTERPOLATION"] = False # Default is to not perform linear interpolation


    ####################################
    #### THIS KEYWORD IS NOT TESTED AND IS NO LONGER NEEEDED ! ####
    try:
        tmp = DYN_PROPERTIES["CHECK_TRIVIAL_CROSSING"]
    except KeyError:
        DYN_PROPERTIES["CHECK_TRIVIAL_CROSSING"] = False # Default is not to check for trivial crossings
    ####################################

    #if ( DYN_PROPERTIES["FUNCTIONAL"] in ["DFTB", "DFTBA"] ):
    #    try:
    #        tmp = DYN_PROPERTIES["BASIS_SET"]
    #    except KeyError:
    #        # For DFTB, we do not use a BASIS_SET definition.
    #        DYN_PROPERTIES["BASIS_SET"] = None
    #    if ( DYN_PROPERTIES["BASIS_SET"] != None ):
    #        # We could also just set to None anyway and continue the calculation. 
    #        #    The current way, however, keeps the user accountable for what they are doing. I like this. ~BMW
    #        print("'BASIS_SET' must not be specified for DFTB or DFTBA Hamiltonians.")
    #        exit()




    print("Input looks good.")















    # Tell user what the code thinks it should run. For double checking the intentions.
    print("Calculation Goal:")
    if ( DYN_PROPERTIES["BOMD"] == True ): # BOMD
        if ( DYN_PROPERTIES["ISTATE"] == 0 ): # BOMD in GS
            print( "\tSQD will perform BOMD in the ground (S0) electronic state." )
            if ( DYN_PROPERTIES['NStates'] >= 2 ):
                print("\nWarning: USER requested multiple electronic states.")
                print("Warning: We will simply be tracking the excited state energies.")
                print("Warning: Set NSTATES = 1 for GS BOMD without calculating excited states.\n")
        if ( DYN_PROPERTIES["ISTATE"] != 0 ): # BOMD in ES
            print( f"\tSQD will perform BOMD in an excited (S{DYN_PROPERTIES['ISTATE']}) electronic state." )
            if ( DYN_PROPERTIES['NStates'] > DYN_PROPERTIES['ISTATE']+1 ):
                print(f"\nWarning: User requested more excited states (NSTATES = {DYN_PROPERTIES['NStates']}) than necessary to perform BOMD in S{DYN_PROPERTIES['ISTATE']}")
                print(f"Warning: Set NSTATES = {DYN_PROPERTIES['ISTATE']+1} for minimal computation time.\n")
    else: # BOMD was not manually set
        if ( DYN_PROPERTIES["CPA"] == False ):
            if ( DYN_PROPERTIES['NStates'] == 1 and DYN_PROPERTIES['ISTATE'] == 0 ): # BOMD in GS -- implicit
                print("User requested 1 electronic state and an initial state of S0 --> BOMD in GS")
                print("\t--> Setting BOMD = True")
                DYN_PROPERTIES["BOMD"] = True
            else: # NAMD
                print( f"\tSQD will perform NAMD in using multiple electronic states")
                print( "\t\tElectronic states in NAMD:", [ f"S{j}" for j in range(DYN_PROPERTIES['NStates']) ] )
                print( f"\t\tInitial Electronic State = S{DYN_PROPERTIES['ISTATE']}.")
        elif ( DYN_PROPERTIES["CPA"] == True ):
            if ( DYN_PROPERTIES['NStates'] == 1 ): # BOMD in GS -- implicit
                print( f"WARNING:\n\tUser chose to do CPA but only one state was given." )
                print( "\tPlease use BOMD keyword for this style of job." )
                exit()
            else:
                print( f"\tSQD will perform NAMD in using multiple electronic states")
                print( "\t\tElectronic states in NAMD:", [ f"S{j}" for j in range(DYN_PROPERTIES['NStates']) ] )
                print( f"\t\tInitial Electronic State = S{DYN_PROPERTIES['ISTATE']}.")
                print( f"\t\tClassical Path Approximation (CPA) for nuclear forces.")


    return DYN_PROPERTIES