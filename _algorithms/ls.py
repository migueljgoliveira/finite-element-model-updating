# IMPORT PACKAGES
from scipy.optimize import least_squares

# IMPORT MODULES

# ------------------------------------------------------------------------ ls()
def ls(methods,ls,opti,options,tests,flag):

    # SET BOUNDS
    if ls.usebounds and ls.method in ['trf','dogbox']:
        bounds = []
        for n in range(0,len(opti.lower)):
            bounds.append([opti.lower[n],opti.upper[n]])
        ls.bounds = bounds

    if ls.x_scale == 'array':
        ls.x_scale = (opti.lower+opti.upper)/2

    result = least_squares(methods,
                           args=(opti,options,tests,flag),
                           x0=opti.initial,
                           jac=ls.jac,
                           bounds=ls.bounds,
                           method=ls.method,
                           ftol=ls.ftol,
                           xtol=ls.xtol,
                           gtol=ls.gtol,
                           x_scale=ls.x_scale,
                           loss=ls.loss,
                           f_scale=ls.f_scale,
                           diff_step=ls.diff_step,
                           tr_solver=ls.tr_solver,
                           tr_options=ls.tr_options,
                           jac_sparsity=ls.jac_sparsity,
                           max_nfev=ls.max_nfev,
                           verbose=ls.verbose)

    return result