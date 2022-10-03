"""
    description
"""

# IMPORT PACKAGES -------------------------------------------------------------
import os
import numpy as np
from shutil import rmtree

# IMPORT MODULES --------------------------------------------------------------
import _output

# -------------------------------------------------- method - OUTPUT DATA FILES
def output_datafiles(opti,options,run):

    # OUTPUT TEXT FILES
    _output.txt(opti,options,run)

    # OUTPUT NUMPY FILES
    _output.npy(opti,options,run)

    # OUTPUT MAT FILES
    _output.mat(opti,options,run)

    return

# --------------------------------------------- main method - OUTPUT SIMULATION
def output_simulation(out,options,tests,run):

    path = 'Output/sim/'+str(run)+'/'
    if run == -1:
        path = 'Output/sim/'
    

    strainfemu = out.fieldoutput.femu.strain
    stressfemu = out.fieldoutput.femu.stress
    forcefemu = out.fieldoutput.femu.force

    for test in tests:
        name = test.info.name
        fpath = path + name + '/'
        try: 
            os.mkdir(fpath)
        except:
            pass
        
        incmax = len(test.data.time)
        with open(fpath+'FEMU_NumStrain_{:s}.dat'.format(name),'a') as f:
            for i in range(0,incmax):
                x = strainfemu[name][i,:,:]
                np.savetxt(f,x,fmt='%+.6e',delimiter=' ')
                f.write('\n')
        
        with open(fpath+'FEMU_NumStress_{:s}.dat'.format(name),'a') as f:
            for i in range(0,incmax):
                x = stressfemu[name][i,:,:]
                np.savetxt(f,x,fmt='%+.6e',delimiter=' ')
                f.write('\n')

        with open(fpath+'FEMU_NumForce_{:s}.dat'.format(name),'a') as f:
            x = forcefemu[name]
            np.savetxt(f,x,fmt='%+.6e',delimiter=' ')
            f.write('\n')

    return

# ------------------------------------------------------------ create_folders()
def create_folders():

    for i in ['txt','npy','mat','sim']:
        try:
            path = f'Output/{i}/'
            os.mkdir(path)
        except:
            pass

# ------------------------------------------------------------ delete_folders()
def delete_folders():

    for i in ['txt','npy','mat','sim']:
        try:
            path = f'Output/{i}/'
            rmtree(path)
        except:
            pass


# -------------------------------------------------------- main method - OUTPUT
def main(out,options,tests,run):

    create_folders()
    if options.info.calibration == False:
        run = -1

    #output_simulation(out,options,tests,run)

    if options.info.calibration == True:
        output_datafiles(out,options,run+1)

    return