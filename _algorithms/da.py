# IMPORT PACKAGES
from scipy.optimize import dual_annealing

# IMPORT MODULES

# ------------------------------------------------------------------------ da()
def da(methods,da,opti,options,tests,flag):

    bounds = []
    for n in range(0,len(opti.lower)):
        bounds.append([opti.lower[n],opti.upper[n]])

    result = dual_annealing(methods, 
                            bounds=bounds,
                            args=(opti,options,tests,flag),
                            maxiter=da.maxiter,
                            local_search_options=da.local_search_options,
                            initial_temp=da.initial_temp,
                            restart_temp_ratio=da.restart_temp_ratio,
                            visit=da.visit,
                            accept=da.accept,
                            maxfun=da.maxfun,
                            no_local_search=da.no_local_search,
                            x0=da.x0)

    return result