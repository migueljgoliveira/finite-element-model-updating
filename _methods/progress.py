# IMPORT PACKAGES -------------------------------------------------------------
import numpy as np
import multiprocessing as mp

def progress(x,objfunc):
    
    # NUMBER OF WORKERS IN PARALLEL PROCESSING
    workers = mp.Process()._identity
    if len(workers) == 2:
        worker = str(mp.Process()._identity[0])
    else:
        worker = '-1'

    if worker == '-1':
        f1 = 'Output/tmp/Var_Evol.dat'
        f2 = 'Output/tmp/Obj_Evol.dat'
    else:
        f1 = 'Output/tmp/Var_Evol_{:s}.dat'.format(worker)
        f2 = 'Output/tmp/Obj_Evol_{:s}.dat'.format(worker)

    with open(f1,'a') as f:
        np.savetxt(f,[x],fmt='%.12e',delimiter=' ',newline='\n')

    with open(f2,'a') as f:
        np.savetxt(f,[objfunc],fmt='%.12e',delimiter=' ',newline='\n')