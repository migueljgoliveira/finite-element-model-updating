"""
    description
"""

# IMPORT PACKAGES -------------------------------------------------------------
import numpy as np
import multiprocessing as mp
from functools import partial
from pprint import pprint as pp

# IMPORT MODULES --------------------------------------------------------------
import _methods

# --------------------------------------------------- main method - METHODS AUX
def methods_aux(test,variables,opti,options,flag):

    fvec = []
    objfunc = []

    weight = test.info.weight
    # FEMU
    if options.method.id == 1:
        femufvec,femuobjfunc = _methods.femu(variables,opti,options,test,flag)

        fvec.append(femufvec*weight)

        objfunc.append(femuobjfunc*weight)

    # VFM
    elif options.method.id == 2:
        vfmfvec,vfmobjfunc = _methods.vfm(variables,opti,options,test,flag)

        fvec.append(vfmfvec*weight)

        objfunc.append(vfmobjfunc*weight)

    # FEMU & VFM
    elif options.method.id == 3:
        femufvec,femuobjfunc = _methods.femu(variables,opti,options,test,flag)
        vfmfvec,vfmobjfunc = _methods.vfm(variables,opti,options,test,flag)

        fvec.append(femufvec*weight)
        fvec.append(vfmfvec*weight)

        objfunc.append(femuobjfunc*weight)
        objfunc.append(vfmobjfunc*weight)

    return fvec,objfunc

# ------------------------------------------------------- main method - METHODS
def main(x,opti,options,tests,flag):

    variables = np.array(x)
    ntests = len(tests)

    # VARIABLES TRANSFORMATION
    if options.info.calibration:
        if options.optimisation.id in [3,4]:
            if options.optimisation.algo.transform:
                lb = opti.lower
                ub = opti.upper
                for n in range(len(variables)):
                    if variables[n] >= 1.0:
                        variables[n] = 1 + (ub[n]-1)*(1-np.exp((1-variables[n])/(ub[n]-1)))
                    else:
                        variables[n] = 1 + (lb[n]-1)*(1-np.exp((1-variables[n])/(lb[n]-1)))

    # DE-NORMALISE VARIABLES
    if options.info.calibration:
        # NORMALISED BY BOUND CONSTRAINTS
        if options.optimisation.id in [1,2]:
            lb = opti.ref['lower']
            ub = opti.ref['upper']
            variables = variables*(ub-lb)+lb
        # NORMALISED BY INITIAL SET
        elif options.optimisation.id in [3,4]:
            variables = variables*opti.ref['initial']

    fvec = []
    objfunc = []

    # EXECUTE TESTS IN PARALLEL
    if options.optimisation.parallel and ntests > 1:
        pool = mp.Pool(ntests)
        part = partial(methods_aux,variables=variables,opti=opti,
                                   options=options,flag=flag)
        res = pool.map(part,tests)

        for i in range(ntests):
            fvec.append(res[i][0])
            objfunc.append(res[i][1])

    # EXECUTE TESTS IN SEQUENCE
    else:
        for test in tests:
            res = methods_aux(test,variables,opti,options,flag)

            fvec.append(res[0])
            objfunc.append(res[1])

    fvec = np.concatenate(fvec,axis=1)
    objfunc = np.concatenate(objfunc,axis=0)

    if options.info.calibration:
        if options.optimisation.result == 'objfunc':
            result = sum(objfunc)
        elif options.optimisation.result == 'fvec':
            result = np.concatenate(fvec)
    else:
        result = sum(objfunc)

    _methods.progress(variables,objfunc)

    return result
