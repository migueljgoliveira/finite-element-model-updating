# IMPORT PACKAGES -------------------------------------------------------------

# IMPORT MODULES --------------------------------------------------------------
import _algorithms

# ---------------------------------------------------- main method - ALGORITHMS
def select(methods,opti,options,tests,flag):

    algo = options.optimisation.algo

    # DIFFERENTIAL EVOLUTION
    if options.optimisation.id == 1:
        result = _algorithms.de(methods,algo,opti,options,tests,flag)

    # DUAL ANNEALING
    elif options.optimisation.id == 2:
        result = _algorithms.da(methods,algo,opti,options,tests,flag)

    # LEAST-SQUARES
    elif options.optimisation.id == 3:
        result = _algorithms.ls(methods,algo,opti,options,tests,flag)
    
    # MINIMIZE
    elif options.optimisation.id == 4:
        result = _algorithms.mini(methods,algo,opti,options,tests,flag)

    return result
