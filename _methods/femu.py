# IMPORT PACKAGES
import os
import time
import signal
import numpy as np
import multiprocessing as mp
from subprocess import Popen,PIPE,TimeoutExpired
from shutil import rmtree,copytree,ignore_patterns

# IMPORT MODULES
import _fortran

# ---------------------------------------------------------------------- femu()
def femu(x,opti,options,test,flag):

    name = test.info.name

    # NUMBER OF WORKERS IN PARALLEL PROCESSING
    worker = '-1'
    if options.optimisation.parallel == False:
        workers = mp.Process()._identity
        if len(workers) == 2:
            worker = str(mp.Process()._identity[0])

    # SET PATH OF ABAQUS EXECUTION
    if worker == '-1':
        path = 'Computation\{:s}'.format(name)
    else:
        path = 'Computation\Worker{:s}\{:s}'.format(worker,name)

    # REMOVE TEST'S OPTIMISATION FOLDER
    try:
        rmtree(path)
    except:
        pass

    # COPY TEST'S INPUT FOLDER
    try:
        copytree('Input\{:s}'.format(name),path,          ignore=ignore_patterns('ExpData.inp','Model.inp'))
    except:
        pass

    # WRITE NEW PARAMETERS TO PARAM.INP
    try:
        with open('{:s}\Param.inp'.format(path),'w') as f:
            f.write('*Parameter')
            for i in range(0,len(x)):
                f.write('\n')
                f.write('V{:d}={:.12f}'.format(i+1,x[i]))
    except:
        pass
    
    # CHANGE DIRECTORY
    cwd = os.getcwd()
    os.chdir(path)

    # EXECUTE ABAQUS JOB
    cpus = str(options.method.femu.cpus)
    umat = options.method.femu.umat
    user = ''
    if umat != None:
        user = f' user={umat}'	
        
    fout = open('abaqus.log','w')

    success = False
    trial = 0
    while not success:
        try:
            cmd = f'abaqus job={name}{user} int cpus={cpus}'
            p = Popen(cmd,cwd=os.getcwd(),stdout=fout,stderr=fout,shell=True)
            p.communicate()
        except TimeoutExpired as e:
            pid = p.pid
            os.kill(pid,signal.SIGTERM)
            os.system('taskkill /F /PID {:d}'.format(pid))
            #cmdkill = 'abaqus job={:s} terminate'.format(name)
            #_ = Popen(cmdkill,cwd=os.getcwd(),stdout=fout,stderr=fout,shell=True).communicate()
            #time.sleep(30)
            #p.kill()
            #p.communicate()
            time.sleep(30)

            for f in os.listdir(os.getcwd()):
                if f.startswith(name) and not f.endswith('inp'):
                    os.remove(f)
            time.sleep(30)

        # EXECUTE ABAQUS PYTHON
        try:
            cmd = 'abaqus python AbaqusPythonScript.py'
            p = Popen(cmd,cwd=os.getcwd(),stdout=fout,stderr=fout,shell=True)
            p.communicate()
        except:
            pass

        # VERIFY NUMERICAL DATA FILES
        if os.path.isfile('NumForce.dat') and os.path.isfile('NumDisplacement.dat'):
            success = True

        trial = trial + 1
        if trial > 2:
            break

    fout.close()
    
    # COPY FEMU (F2PY) ARGUMENTS (I)
    incmax = len(test.data.time)
    nelems = len(test.model.elements)
    ntests = len(options.info.name)
    norm = options.method.femu.norm
    residual = options.method.residual
    femuresidual = options.method.femu.residual
    exptime = test.data.time
    expforce = test.data.force
    expstrain = test.data.strain
    maxexpforce = np.max(test.data.force)
    maxexpstrain = np.max(test.data.strain)

    # COMPUTE NUMBER OF RESIDUALS
    if femuresidual == 0:
        m = incmax
    elif femuresidual == 1:
        m = incmax + incmax
    elif femuresidual == 2:
        m = (incmax*nelems*3)+incmax

    if success:

        # READ NUMERICAL DATA
        auxnumtime,auxnumforce = np.loadtxt('NumForce.dat',delimiter=',',unpack=True)
        incnum = len(auxnumtime)
        nnodes = len(test.model.coordinates)
        auxnumdispl = np.zeros((incnum,nnodes,2))

        for i in range(incnum):
            auxnumdispl[i,:,:] = np.loadtxt('NumDisplacement.dat',  delimiter=',',skiprows=i*nnodes,max_rows=nnodes)

        # COPY DISPLACEMENTS_TO_STRAINS (F2PY) ARGUMENTS
        nelems = len(test.model.elements)
        nodpelem = len(test.model.elements[0])
        incnum = len(auxnumtime)
        nnodes = len(test.model.coordinates)
        elems = test.model.elements
        nodes = test.model.coordinates
        displ = auxnumdispl

        # CALL DISPLACEMENTS_TO_STRAINS (F2PY)
        auxnumstrain,_,_,_,_,_ = _fortran.displacements_to_strains(nelems,nodpelem,incnum,nnodes,elems,nodes,displ,False)

        # COPY FEMU (F2PY) ARGUMENTS (II)
        incnum = len(auxnumtime)
        numtime = auxnumtime
        auxnumforce = auxnumforce
        auxnumstrain = auxnumstrain
        
        # CALL FEMU (F2PY)
        fvec,objfunc = _fortran.femu(m,incmax,incnum,nelems,ntests,norm,residual,femuresidual,exptime,expforce,expstrain,numtime,auxnumforce,auxnumstrain,maxexpforce,maxexpstrain)
    else:

        fvec = np.ones(m)
        objfunc = 1.0

    # CHANGE TO OLD DIRECTORY
    os.chdir(cwd)

    # NORMALIZE BY FEMU WEIGHT
    weight = options.method.femu.weight

    fvec = fvec*weight
    objfunc = objfunc*weight

    return fvec,objfunc