#!/usr/bin/python
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import re
import shutil
import json
import time


def dump_data(filename, obj):
    """
    dump data in to disk in a special format.
    wrap internal storage method.
    """
    dump_json(filename, obj)        
    
    return
    
def load_data(filename):
    """
    load data into memory from disk
    wrap internal storage method.
    """
    if os.path.isfile(filename):
        obj = load_json(filename)
    else:
        print ("CANNOT FIND file: %s" % filename)
        sys.exit(1)
    
    return obj

def dump_json(filename, obj, encode='utf-8'):
    """
    dump an object in json format.
    """
    json.encoder.FLOAT_REPR = lambda f: format("%.18g" % f)
    fp = open(filename, mode='w')
    my_str = json.dumps(obj, encoding=encode, indent=2)
    fp.write(my_str)
    fp.close()
    
    return

def load_json(filename, encode='utf-8'): 
    """
    load an object in json format.
    """        
    json.encoder.FLOAT_REPR = lambda f: format("%.18g" % f)
    fp = open(filename, mode='r')
    obj = json.load(fp, encoding='utf-8')
    
    return obj     

class gau_nac:
    """
    calc. nac data.
    gaussian.chk
    """
    def __init__(self, DYN_PROPERTIES,  config = {}):
        """
        very ok nac
        """
        self.DYN_PROPERTIES = DYN_PROPERTIES
        
        self.directory = {'work': "TD_NEW_S1", \
                          'work_prev': "TD_OLD_S1", \
                          "DIMER": "DIMER", \
                          "OVERLAP": "OVERLAP"  \
                          }
        self.files = {'dimension': "dimension.json"}
        
        if config != {}:
            root_dir = config['root']
            dirs = config['dirs']
            files = config['files']            
            
            # working directory & files >>>
            self.directory = {}
            self.directory['root'] = root_dir
            self.directory['home'] = root_dir + "/" + dirs['home']
            
            self.directory['work'] = self.directory['home'] + "/" + dirs['work']
            self.directory['work_prev'] = self.directory['home'] + "/" + dirs['work_prev']
            self.directory["DIMER"] = self.directory['home'] + "/" + dirs["DIMER"]
            self.directory["OVERLAP"] = self.directory['home'] + "/" + dirs["OVERLAP"]
            self.files = {}
            self.files["dimension"] = files['dimension']
            
            # run the job directly
            self.worker()
                        
        return
        

    def prepare(self):
        """
        first, prepare work dir; then, the necessary files.
        """
        work_dir = self.directory["OVERLAP"]
        if os.path.exists(work_dir):
	        shutil.rmtree(work_dir)
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
   
        sourceFile = self.directory['work_prev'] + '/mo.dat'
        destFile = self.directory["OVERLAP"] + '/mo_1.dat'
        shutil.copy2(sourceFile, destFile)
 
        sourceFile = self.directory['work_prev'] + '/ci.dat'
        destFile   = self.directory["OVERLAP"] + '/ci_1.dat'
        shutil.copy2(sourceFile, destFile)
 
        sourceFile = self.directory['work'] + '/mo.dat'
        destFile =   self.directory["OVERLAP"] + '/mo_2.dat'
        shutil.copy2(sourceFile, destFile)

        sourceFile = self.directory['work'] + '/ci.dat'
        destFile =   self.directory["OVERLAP"] + '/ci_2.dat'
        shutil.copy2(sourceFile, destFile)

        #sourceFile = self.directory['work'] + '/qm_results.dat'
        #destFile   = self.directory["OVERLAP"]  + '/qm_results.dat'
        #shutil.copy2(sourceFile, destFile)
        
        sourceFile = self.directory['work'] + '/' + self.files['dimension']
        destFile   = self.directory["OVERLAP"]  + '/' + self.files['dimension']
        shutil.copy2(sourceFile, destFile)
        
        sourceFile = self.directory["DIMER"] + '/ao_overlap.dat'
        destFile =   self.directory["OVERLAP"]  + '/ao_overlap.dat'
        shutil.copy2(sourceFile, destFile)

        os.chdir(work_dir)
        
        # load internal data.
        filename = self.files['dimension']
        dim = load_data(filename)
        
        n_atom = dim['n_atom']      # Number of atom
        n_state = dim['n_state']    # Number of states
        n_ao = dim['n_basis']       # number of basis functions
        n_occ = dim['nocc_allA']    # number of occupied orbitals        
        tmp = len( open('ci_1.dat','r').readlines()[1:] )
        n_csf = tmp // (n_state-1) # number of excited slater determinents per excited state       

        fileout1=open('main_overlap_slater_input','w')
        fileout1.write('                        read (*,*)  \n')
        fileout1.write(''+str(n_atom)+'               read (*,*) n_atom \n')
        fileout1.write(''+str(n_ao)+'                 read (*,*) n_ao \n')
        fileout1.write(''+str(n_occ)+'               read (*,*) n_ele_alpha \n')
        fileout1.write(''+str(n_occ)+'               read (*,*) n_ele_beta \n')
        fileout1.write('                        read (*,*)  \n')
        fileout1.write(''+str(n_state)+'               read (*,*) n_state \n')
        fileout1.write(''+str(n_csf)+'               read (*,*) n_csf \n')
        fileout1.write('                        read (*,*)  \n')
        fileout1.write('1                       read (*,*)  type_input  \n')
        fileout1.write('ci_1.dat                read (*,*)  filename_input1  \n')
        fileout1.write('ci_2.dat                read (*,*)  filename_input2  \n')
        fileout1.write('overlap.dat             read (*,*)  filename_input2  \n')
        fileout1.write('                        read (*,*)  \n')
        fileout1.write('0                       read (*,*) output_level  \n')
        fileout1.write('ci_overlap.dat          read (*,*) filename_output  \n')
        fileout1.close()

        #print("OVERLAP Extracted.")
        
        return
###

    def run(self):
        """
        call another standalone program to deal with nac.
        """
        #os.system("main_overlap_slater.exe  < main_overlap_slater_input")
        os.system(f"{self.DYN_PROPERTIES['FORTRAN_CI_CODE']}  < main_overlap_slater_input")
        return
        
    def dump(self):
        """
            dump necessary data of nac
        """
        # open file & read data.
        fp = open('qm_results.dat','r')
        qms = fp.readlines()
        fp.close()
        
        # num. of states
        for cur_line in qms:
            i_find_n_state = re.search ('Number of electronic states', cur_line)
            if i_find_n_state is not None:
                n_state = int (cur_line.split()[-1]) 
                               
        # open & read wf-overlap
        fp = open('wavefuction_overlap.dat','r') 
        ci_overlap = fp.readlines()        
        fp.close()         
        
        # dump data.
        n_line = len(qms)
        fp = open('qm_result_update.dat','w')
        for i_line in range(n_line - n_state*n_state - 2): # n_state^2 + 2 is wf data.
            fp.write(''+str(qms[i_line][:-1]) + '  \n')
            
        fp.write('Wave-function overlap between R and R+dR \n')          
        i_line=1 
        
        for i_state in range(n_state):
            for j_state in range(n_state):    
                fp.write('S'+str(i_state)+'   S'+str(j_state)+'    '+ str(ci_overlap[i_line].split()[-1])+'  \n')
 
                i_line = i_line + 1
        fp.write('----------------------------------------------')
        fp.close()           

        return

    #def finilize(self):
        """
        finish the current step & prepare for the following step
        """
        # file man. & back up
        #shutil.copyfile("./qm_result_update.dat", "./qm_results.dat")
        #shutil.copyfile("./qm_results.dat", "../../qm_results.dat")          
        #   Go back to directory of dynamics work
        #os.chdir("../")     
           
        return

    def dump_braden(self):

        # open & read wf-overlap
        f = open('wavefuction_overlap.dat','r') 
        lines = f.readlines()[1:]
        f.close()
        NStates = int(lines[-1].split()[0])
        CI_overlap = np.zeros(( NStates, NStates ))
        for line in lines:
            t = line.split()
            j = int( t[0] ) - 1
            k = int( t[1] ) - 1
            CI_overlap[j,k] = float( t[2] )
        np.savetxt("wavefunction_overlap_MAT.dat", CI_overlap)
        plt.imshow( np.abs(CI_overlap), origin='lower', cmap="hot_r" )
        plt.colorbar(pad=0.01)
        plt.xlabel("Electronic State Index",fontsize=15)
        plt.ylabel("Electronic State Index",fontsize=15)
        plt.tight_layout()
        plt.savefig("wavefunction_overlap_MAT.jpg",dpi=600)
        plt.clf()
 
        # Recall, we compute one additional state in TD-DFT. Do not save it here.
        if ( self.DYN_PROPERTIES["MD_STEP"] >= 2 ):
            self.DYN_PROPERTIES["OVERLAP_OLD"] = (self.DYN_PROPERTIES["OVERLAP_NEW"])
        self.DYN_PROPERTIES["OVERLAP_NEW"] = (CI_overlap[:,:])[:-1,:-1]




    def worker(self):
        """
        prepare; run; dump; finilize
        """
        if ( self.DYN_PROPERTIES["MD_STEP"] == 1 ):
            self.DYN_PROPERTIES["FORTRAN_CI_CODE"] = f"{self.DYN_PROPERTIES['SQD_HOME_PATH']}/src/WFN_OVERLAP/FORTRAN/main_overlap_slater.exe"


        T0 = time.time()
        self.prepare()
        print("gau_nac prepare", round(time.time() - T0,2), "s")
        T0 = time.time()
        self.run()
        print("WFN OVERLAP FORTRAN TOOK", round(time.time() - T0,2), "s")
        #self.dump()
        T0 = time.time()
        self.dump_braden()
        print("gau_nac dump", round(time.time() - T0,2), "s")
        os.chdir("../")
        
        return self.DYN_PROPERTIES

# main program.
if __name__ == "__main__":  
    DYN_PROPERTIES = {"test":None}  
    n = gau_nac(DYN_PROPERTIES) 
    DYN_PROPERTIES = n.worker()

   
     
     
      
