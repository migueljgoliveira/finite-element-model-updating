# IMPORT PACKAGES -------------------------------------------------------------
import os
import shutil
import numpy as np
from scipy._lib._util import check_random_state

# IMPORT MODULES --------------------------------------------------------------
import _output
import _methods
import _algorithms

# -------------------------------------------------------- class - OPTIMISATION
class Optimisation:
    def __init__(self,options):
        self.ref = options.material.parameters
        self.initial = self.ref['initial']
        self.lower = self.ref['lower']
        self.upper = self.ref['upper']

        self.objfunc = ObjectiveFunction()
        self.variables = Variables()

        self.best = Best()

        self.fieldoutput = FieldOutput()

    def normalised(self,options):
        # NORMALISED BY BOUND CONSTRAINTS
        if options.optimisation.id in [1,2]:
            lb = self.lower
            ub = self.upper
            self.initial = (self.ref['initial']-lb)/(ub-lb)
            self.lower = np.zeros(len(lb))
            self.upper = np.ones(len(ub))
        # NORMALISED BY INITIAL SET
        elif options.optimisation.id in [3,4]:
            self.initial = self.initial/self.ref['initial']
            self.lower = self.lower/self.ref['initial']
            self.upper = self.upper/self.ref['initial']

class ObjectiveFunction:
    def __init__(self):
        self.best = Method()
        self.all = Method()
    
class Variables:
    def __init__(self):
        self.best = None
        self.all = None

class Method:
    def __init__(self):
        self.vfm = None
        self.femu = None
        self.total = None

class Best:
    def __init__(self):
        self.evolution = None
        self.evaluation = None
        self.variables = None
        self.objfunc = Method()

class FieldOutput:
    def __init__(self):
        self.femu = FEMU()
        self.vfm = VFM()

class FEMU:
    def __init__(self):
        self.strain = dict()
        self.force = dict()
        self.stress = dict()
        self.displacement = None

class VFM:
    def __init__(self):
        self.stress = dict()

# -------------------------------------------- main method - LOAD FIELD OUTPUTS
def load_field_output(opti,options,tests):

    path = 'Output/tmp/'
    for test in tests:
        name = test.info.name
        incmax = len(test.data.time)
        nelems = len(test.model.elements)
        if options.method.id in [1,3]:
            opti.fieldoutput.femu.strain[name] = np.zeros((incmax,nelems,3))
            opti.fieldoutput.femu.stress[name] = np.zeros((incmax,nelems,3))
            opti.fieldoutput.femu.force[name] = np.zeros(incmax)
            for i in range(0,incmax):
                data = np.loadtxt(path+'FEMU_NumStrain_{:s}.dat'.format(name),skiprows=i*nelems,max_rows=nelems)
                opti.fieldoutput.femu.strain[name][i,:,:] = data
                data = np.loadtxt(path+'FEMU_NumForce_{:s}.dat'.format(name))
                opti.fieldoutput.femu.force[name][:] = data
                data = np.loadtxt(path+'FEMU_NumStress_{:s}.dat'.format(name),skiprows=i*nelems,max_rows=nelems)
                opti.fieldoutput.femu.stress[name][i,:,:] = data

        if options.method.id in [2,3]:
            opti.fieldoutput.vfm.stress[name] = np.zeros((incmax,nelems,3))
            for i in range(0,incmax):
                data = np.loadtxt(path+'VFM_NumStress_{:s}.dat'.format(name),skiprows=i*nelems,max_rows=nelems)
                opti.fieldoutput.vfm.stress[name][i,:,:] = data

    return opti

# ----------------------------------------------------- main method - FIND BEST
def find_best(opti):
    objfunc = opti.objfunc.all.total
    check = objfunc[0]
    if not isinstance(check,float):
        objfunc = np.sum(objfunc,axis=1)

    idxs = [0]
    current_best = objfunc[0]
    for n in range(0,len(objfunc)):
        if objfunc[n] < current_best:
            current_best = objfunc[n]
            idxs.append(n)
    idxs.append(len(objfunc)-1)

    opti.best.evolution = np.zeros(len(opti.variables.all),dtype=int)
    opti.variables.best = np.full_like(opti.variables.all, 0.0)
    opti.objfunc.best.total = np.full_like(opti.objfunc.all.total, 0.0)
    opti.objfunc.best.vfm = np.full_like(opti.objfunc.all.vfm, 0.0)
    opti.objfunc.best.femu = np.full_like(opti.objfunc.all.femu, 0.0)

    for n in range(0,len(idxs)-1):
        i = idxs[n]
        f = idxs[n+1]+1

        opti.best.evolution[i:f] = n+1
        opti.variables.best[i:f] = opti.variables.all[i]

        opti.objfunc.best.femu[i:f] = opti.objfunc.all.femu[i]
        opti.objfunc.best.vfm[i:f] = opti.objfunc.all.vfm[i]
        opti.objfunc.best.total[i:f] = opti.objfunc.all.total[i]

    opti.best.evaluation = idxs[-2]
    opti.best.variables = opti.variables.best[-1]
    opti.best.objfunc.vfm = opti.objfunc.best.vfm[-1]
    opti.best.objfunc.femu = opti.objfunc.best.femu[-1]
    opti.best.objfunc.total = opti.objfunc.best.total[-1]

    return

# ----------------------------------------------------- main method - LOAD DATA
def load_data_opti():

    # List files used in parallel processing
    files1 = []
    files2 = []
    for f in os.listdir('Output/tmp/'):
        if f.startswith('Var_Evol_'):
            files1.append(f)
        elif f.startswith('Obj_Evol_'):
            files2.append(f)

    files1.sort()
    files2.sort()

    # Load data for single processing
    if len(files1) == 0:
        with open('Output/tmp/Var_Evol.dat','r') as f:
            variables = np.loadtxt(f,delimiter=' ',ndmin=2)

        with open('Output/tmp/Obj_Evol.dat','r') as file:
            objfunc = np.loadtxt(file,delimiter=' ',ndmin=2)

    # Load data for multiprocessing
    else:
        auxvariables = []
        for file1 in files1:
            with open('Output/tmp/{:s}'.format(file1),'r') as f:
                auxvariables.append(np.loadtxt(f,delimiter=' '))

        variables = []
        nmax = max(map(len,auxvariables)) 
        n = 0
        while n < nmax:
            for i in range(0,len(auxvariables)):
                try:
                    variables.append(auxvariables[i][n])
                except:
                    pass
            n += 1

        try:
            with open('Output/tmp/Var_Evol.dat','r') as f:
                auxvariables = np.loadtxt(f,delimiter=' ')
            if auxvariables.ndim in [0,1]:
                variables.append(auxvariables)
            else:
                for item in auxvariables:
                    variables.append(item)
        except:
            pass

        variables = np.array(variables)

        auxobjfunc = []
        for file2 in files2:
            with open('Output/tmp/{:s}'.format(file2),'r') as f:
                auxobjfunc.append(np.loadtxt(f,delimiter=' ',ndmin=2))

        objfunc = []
        nmax = max(map(len,auxobjfunc)) 
        n = 0
        while n < nmax:
            for i in range(0,len(auxobjfunc)):
                try:
                    objfunc.append(auxobjfunc[i][n])
                except:
                    pass
            n +=1

        try:
            with open('Output/tmp/Obj_Evol.dat','r') as f:
                auxobjfunc = np.loadtxt(f,delimiter=' ',ndmin=1)
            if auxobjfunc.ndim in [1]:
                objfunc.append(auxobjfunc)
            else:
                for item in auxobjfunc:
                    objfunc.append(item)
        except:
            pass
        
        objfunc = np.array(objfunc)

    return variables,objfunc

# ------------------------------------------------ main method - STORE PROGRESS
def store_progress(opti,options,tests):

    # Load field output results of best solution
    # opti = load_field_output(opti,options,tests)

    # Load optimisation results from .dat files
    variables,objfunc  = load_data_opti()

    opti.variables.all = variables

    if options.method.id == 1:
        opti.objfunc.all.femu = objfunc
        opti.objfunc.all.vfm = np.full_like(objfunc,np.nan)
        opti.objfunc.all.total = objfunc
    elif options.method.id == 2:
        opti.objfunc.all.femu = np.full_like(objfunc,np.nan)
        opti.objfunc.all.vfm = objfunc
        opti.objfunc.all.total = objfunc
    elif options.method.id == 3:
        ntests = len(options.info.name)
        nevals = len(objfunc[:,0])
        opti.objfunc.all.femu = np.zeros((nevals,ntests))
        opti.objfunc.all.vfm = np.zeros((nevals,ntests))
        for n in range(0,ntests):
            opti.objfunc.all.femu[:,n] = objfunc[:,2*n]
            opti.objfunc.all.vfm[:,n] = objfunc[:,2*n+1]
            opti.objfunc.all.total = objfunc[:,2*n] + objfunc[:,2*n+1]	

    # Find updated current best from all evaluations
    find_best(opti)
    
    return opti

# ------------------------------------------------------ method - PRINT RESULTS
def print_results(result,options,run,runs,log):

    pvars = ['%.12f'%i for i in result.x]
    pvars = ' , '.join(pvars)
    if options.optimisation.id in [1,2,4]:
        cost = result.fun
    else:
        cost = result.cost

    if run == 0:
        log.info(' '.ljust(21)+'> OUTPUT')

    log.info('\n')
    if runs != 1:
        log.info(' '.ljust(10)+'        RUN: {}'.format(run+1))
    log.info(' '.ljust(10)+'    SUCCESS: {}'.format(result.success))
    log.info(' '.ljust(10)+'    MESSAGE: {:s}'.format(result.message))
    log.info(' '.ljust(10)+'       COST: {:12e}'.format(cost))
    log.info(' '.ljust(10)+'  VARIABLES: {}'.format(pvars))
    try:
        log.info(' '.ljust(10)+' ITERATIONS: {:d}'.format(result.nit))
    except:
        pass
    log.info(' '.ljust(10)+'EVALUATIONS: {:d}'.format(result.nfev))

# -------------------------------------------------- method - OPTIMISATION INIT
def optimisation_init(options,tests,runs,run,log):

    opti = Optimisation(options)

    # NORMALISE VARIABLES
    opti.normalised(options)

    # CREATE OPTIMISATION FOLDER
    try:
        os.mkdir('Computation/')
    except:
        pass

    # REMOVE OLD FOLDERS
    try:
        for f in os.listdir('Computation/'):
            shutil.rmtree('Computation/'+f)
    except:
        pass
    
    # CALL ALGORITHMS
    flag = False
    result = _algorithms.select(_methods.main,opti,options,tests,flag)

    # COMPUTE OUTPUT FIELDS
    flag = True
    _ = _methods.main(result.x,opti,options,tests,flag)

    # STORE RESULTS
    opti = store_progress(opti,options,tests)

    # PRINT RESULTS
    print_results(result,options,run,runs,log)
    
    return opti

# -------------------------------------------------------- method - MULTI-START
def multi_start(runs,options):

    nvars = len(options.material.parameters['initial'])
    lb = options.material.parameters['lower']
    ub = options.material.parameters['upper']

    # RANDOM SAMPLING
    if options.optimisation.algo.multi == 'Random':
        sets = np.random.random((runs,nvars))
    
    # LATIN HYPER CUBE SAMPLING
    elif options.optimisation.algo.multi == 'LHS':
        rng = check_random_state(1)
        segsize = 1.0 / runs
        samples = segsize * rng.uniform(size=(runs,nvars))
        samples += np.linspace(0,1,runs,endpoint=0)[:,np.newaxis]
        sets = np.zeros_like(samples)

        for j in range(nvars):
            order = rng.permutation(range(runs))
            sets[:, j] = samples[order, j]

    # elif options.optimisation.algo.multi == 'MC':
        # sets = lhsmdu.createRandomStandardUniformMatrix(nvars,runs)

    for i in range(len(sets)):
        sets[i] = lb + sets[i]*(ub-lb)

    return sets

# -------------------------------------------------------- method - TEMP FOLDER
def temp_folder(flag):

    path = 'Output/tmp/'

    # CREATE FOLDER
    if flag == 1:
        try:
            shutil.rmtree(path)
        except:
            pass
        try:
            os.mkdir(path)
        except:
            pass

    # DELETE OLD FILES
    elif flag == 2:
        try:
            for f in os.listdir(path):
                os.remove(path+f)
        except:
            pass

    # DELETE FOLDER
    elif flag == 3:
        try:
            shutil.rmtree(path)
        except:
            pass

# -------------------------------------------------- main method - OPTIMISATION
def main(options,tests,log):

    # DELETE OUTPUT FOLDERS
    _output.delete_folders()

    # CREATE TEMPORARY FOLDER
    temp_folder(1)

    # GENERATE MULTI-START SETS
    runs = 1
    if options.optimisation.id in [3,4]:
        if options.optimisation.algo.multistart:
            runs = options.optimisation.algo.nsamples
            sets = multi_start(runs,options)

    # START OPTIMISATION RUNS
    for run in range(runs):
        if options.optimisation.id in [3,4]:
            if options.optimisation.algo.multistart:
                options.material.parameters['initial'] = sets[run,:]

        out = optimisation_init(options,tests,runs,run,log)

        # SAVE OUTPUT
        _output.main(out,options,tests,run)

        # DELETE OLD FILES FROM TEMPORARY FOLDER
        temp_folder(2)

    # DELETE TEMPORARY FOLDER
    temp_folder(3)