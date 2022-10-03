# IMPORT PACKAGES -------------------------------------------------------------
import os
import numpy as np

# IMPORT MODULES --------------------------------------------------------------
import _methods
import _algorithms

# ---------------------------------------------------------- class - SIMULATION
class Simulation:
    def __init__(self):
        self.objfunc = Method()
        self.variables = None

        self.fieldoutput = FieldOutput()

class Method:
    def __init__(self):
        self.vfm = None
        self.femu = None
        self.total = None

class FieldOutput:
    def __init__(self):
        self.femu = FEMU()
        self.vfm = VFM()

class FEMU:
    def __init__(self):
        self.strain = dict()
        self.force = dict()
        self.stress = dict()
        self.displacement = dict()

class VFM:
    def __init__(self):
        self.stress = dict()

# ------------------------------------------------ main method - STORE PROGRESS
def store_progress(sim,options,tests):

    sim = _algorithms.load_field_output(sim,options,tests)	

    # Load optimisation results from .dat files
    with open('Output/tmp/Var_Evol.dat','r') as f:
        variables = np.loadtxt(f,delimiter=' ')

    with open('Output/tmp/Obj_Evol.dat','r') as file:
        objfunc = np.loadtxt(file,delimiter=' ')

    sim.variables = variables

    if options.method.id == 1:
        sim.objfunc.femu = objfunc
        sim.objfunc.vfm = np.full_like(objfunc,np.nan)
        sim.objfunc.total = objfunc
    elif options.method.id == 2:
        sim.objfunc.femu = np.full_like(objfunc,np.nan)
        sim.objfunc.vfm = objfunc
        sim.objfunc.total = objfunc
    elif options.method.id == 3:
        ntests = len(options.info.name)
        sim.objfunc.femu = np.zeros(ntests)
        sim.objfunc.vfm = np.zeros(ntests)
        for n in range(0,ntests):
            sim.objfunc.femu[n] = objfunc[2*n]
            sim.objfunc.vfm[n] = objfunc[2*n+1]
            sim.objfunc.total = objfunc[2*n] + objfunc[2*n+1]

    try:
        for f in os.listdir('Output/tmp/'):
            os.remove('Output/tmp/'+f)
        os.rmdir('Output/tmp/')
    except:
        pass
    
    return sim

# ---------------------------------------------------- main method - SIMULATION
def main(options,tests):

    x = options.material.parameters['initial']

    sim = Simulation()

    # REMOVE TEMPORARY FILES USED
    try:
        for f in os.listdir('Output/tmp/'):
            os.remove('Output/tmp/'+f)
    except:
        pass

    try:
        os.mkdir('Output/tmp/')
    except:
        pass
    
    # CALL METHODS
    flag = True
    _ = _methods.main(x,sim,options,tests,flag)

    # STORE RESULTS
    sim = store_progress(sim,options,tests)

    pvars = ['%.6f'%i for i in sim.variables]
    pvars = ' , '.join(pvars)
    cost = np.sum(sim.objfunc.total)

    print(' '.ljust(21)+'> OUTPUT')
    print('\n')
    print(' '.ljust(10)+'       COST: {:6e}'.format(cost))
    print(' '.ljust(10)+'  VARIABLES: {}'.format(pvars))

    return sim
