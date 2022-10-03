# IMPORT PACKAGES
from scipy.optimize import minimize

# IMPORT MODULES

# ---------------------------------------------------------------------- mini()
def mini(methods,mini,opti,options,tests,flag):

    # SET BOUNDS
    meths = ['Powell','L-BFGS-B','TNC','SLSQP','trust-constr']

    if mini.usebounds and mini.method in meths:
        bounds = []
        for n in range(0,len(opti.lower)):
            bounds.append([opti.lower[n],opti.upper[n]])
        mini.bounds = bounds

    result = minimize(methods,
                      x0=opti.initial,
                      args=(opti,options,tests,flag),
                      method=mini.method,
                      jac=mini.jac,
                      hess=mini.hess,
                      hessp=mini.hessp,
                      bounds=mini.bounds,
                      tol=mini.tol,
                      options=mini.options)

    return result