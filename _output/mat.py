# IMPORT PACKAGES
import numpy as np
from scipy.io import savemat
from os import mkdir,listdir,remove

# IMPORT MODULES

# ------------------------------------------------------------------ mat_save()
def mat_save(opti,options,run):
    path = 'Output/mat/'+str(run)+'/'
    try:
        for f in listdir(path):
            remove(path+f)
    except:
        pass
    
    try:
        mkdir(path)
    except:
        pass
    
    nevals = len(opti.variables.all)
    evals = np.arange(1,nevals+1)
    best = opti.best.evolution

    data = {'variables'  : {
                            'all' : opti.variables.all,
                            'best': opti.variables.best
                           },
            'objfunc'	 : {
                            'all' : {
                                     'total' : opti.objfunc.all.total,
                                     'femu'  : opti.objfunc.all.femu,
                                     'vfm'   : opti.objfunc.all.vfm
                                    },
                            'best': {
                                     'total' : opti.objfunc.best.total,
                                     'femu'  : opti.objfunc.best.femu,
                                     'vfm'   : opti.objfunc.best.vfm
                                    }
                           },
            'evaluations': {
                            'all'  : evals,
                            'best' : best
                           }
           
           }
    savemat(path+'Evol.mat',data,appendmat=True,oned_as='column')