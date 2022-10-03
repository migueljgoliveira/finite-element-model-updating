# IMPORT PACKAGES
from scipy.optimize import Bounds
from scipy.optimize import differential_evolution

# IMPORT MODULES

# ------------------------------------------------------------------------ de()
def de(methods,de,opti,options,tests,flag):

    lb = opti.lower
    ub = opti.upper
    result = differential_evolution(methods,
                                    args=(opti,options,tests,flag),
                                    bounds=Bounds(lb,ub,keep_feasible=True),
                                    strategy=de.strategy,
                                    maxiter=de.maxiter,
                                    popsize=de.popsize,
                                    tol=de.tol,
                                    mutation=de.mutation,
                                    recombination=de.recombination,
                                    disp=de.disp,
                                    polish=de.polish,
                                    init=de.init,
                                    atol=de.atol,
                                    updating=de.updating,
                                    workers=de.workers)

    return result