# IMPORT PACKAGES -------------------------------------------------------------
import numpy as np
from pprint import pprint
from scipy.io import savemat

# IMPORT MODULES --------------------------------------------------------------
import _input
import _fortran

# ---------------------------------------------------------------- method - LOG
def input_log(options,log):

    method = options.method.description
    algorithm = options.optimisation.description

    log.info(' '.ljust(21)+'> INPUT')
    log.info('\n')
    log.info(' '.ljust(10)+'     METHOD: {}'.format(method))
    log.info(' '.ljust(10)+'  ALGORITHM: {}'.format(algorithm))
    log.info('\n\n')

# ------------------------------------------------------------ method - STRAINS
def input_strains(test):
    
    nelems = len(test.model.elements)
    nodpelem = len(test.model.elements[0])
    incmax = len(test.data.time)
    nnodes = len(test.model.coordinates)
    elems = test.model.elements
    nodes = test.model.coordinates
    displ = test.data.displacements

    strain,dfgrd,indfgrd,detdfgrd,rot,area = _fortran.displacements_to_strains(nelems,nodpelem,incmax,nnodes,elems,nodes,displ,True)

    test.data.allocate(strain,dfgrd,indfgrd,detdfgrd,rot,area)

# --------------------------------------------------------- main method - INPUT
def main(logg):

    # load program settings
    options = _input.Options()

    # load tests data
    tests = []
    n = 0
    for test in options.info.name:
        weight = options.info.weight[n]
        tests.append(_input.Test(test,weight))
        n += 1

    # convert displacements to strains
    for test in tests:
        input_strains(test)

    # print info to command line and log file
    input_log(options,logg)

    return options,tests