import numpy as np


def momentOfInertia(DYN_PROPERTIES):
    MOI     = np.zeros((3,3))
    one_mat = np.identity(3)

    XYZ = DYN_PROPERTIES["Atom_coords_new"]

    for at in range(DYN_PROPERTIES["NAtoms"]):
        MOI += DYN_PROPERTIES["MASSES"][at] * \
            ( np.dot(XYZ[at,:],XYZ[at,:])*one_mat - \
              np.outer(XYZ[at],XYZ[at]))
    
    DYN_PROPERTIES["INERTIA"] = MOI
    return DYN_PROPERTIES

def myCrossProduct(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],\
                        a[2]*b[0]-a[0]*b[2],\
                        a[0]*b[1]-a[1]*b[0]])

def angularMomentum(DYN_PROPERTIES):
    L = np.array([0.0,0.0,0.0])
    XYZ = DYN_PROPERTIES["Atom_coords_new"]
    VEL = DYN_PROPERTIES["Atom_velocs_new"]
    for at in range(DYN_PROPERTIES["NAtoms"]):
        L += DYN_PROPERTIES["MASSES"][at] * \
                myCrossProduct(XYZ[at],VEL[at])
    DYN_PROPERTIES["ANG_MOM"] = L
    return DYN_PROPERTIES

def getAngularVelocity(DYN_PROPERTIES):
    DYN_PROPERTIES = angularMomentum(DYN_PROPERTIES)
    DYN_PROPERTIES = momentOfInertia(DYN_PROPERTIES)
    omega = np.dot(np.linalg.inv(DYN_PROPERTIES["INERTIA"]),DYN_PROPERTIES["ANG_MOM"])
    DYN_PROPERTIES["ANG_VEL"] = omega
    return DYN_PROPERTIES

def remove_rotations(DYN_PROPERTIES):
    DYN_PROPERTIES = getAngularVelocity(DYN_PROPERTIES)
    for at in range(DYN_PROPERTIES["NAtoms"]):
        DYN_PROPERTIES["Atom_velocs_new"][at,:] -= \
            myCrossProduct(DYN_PROPERTIES["ANG_VEL"][:],DYN_PROPERTIES["Atom_coords_new"][at,:])
    return DYN_PROPERTIES


def shift_COM(DYN_PROPERTIES):


    COM = np.zeros((3))
    for at in range(DYN_PROPERTIES["NAtoms"]):
        for d in range(3):
            COM[d] += DYN_PROPERTIES["MASSES"][at] * DYN_PROPERTIES["Atom_coords_new"][at,d]

    M_Total = np.sum( DYN_PROPERTIES["MASSES"] )

    DYN_PROPERTIES["COM"] = COM / M_Total
    DYN_PROPERTIES["Atom_coords_new"] -= DYN_PROPERTIES["COM"]

    return DYN_PROPERTIES