# IMPORT PACKAGES
import numpy as np

# IMPORT MODULES

# ---------------------------------------------------------------- Optimisation
class Optimisation:
    def __init__(self):
        self.id = 3
        self.parallel = True

        if self.id == 1:
            self.name = 'de'
            self.description = 'Differential Evolution'
            self.algo = DE()
            self.result = 'objfunc'
        elif self.id == 2:
            self.name = 'da'
            self.description = 'Dual Annealing'
            self.algo = DA()
            self.result = 'objfunc'
        elif self.id == 3:
            self.name = 'ls'
            self.algo = LS()
            self.description = 'Least-Squares - {:s}'.format(self.algo.method)
            self.result = 'fvec'
        elif self.id == 4:
            self.name = 'mini'
            self.algo = MINI()
            self.description = 'Minimize - {:s}'.format(self.algo.method)
            self.result = 'objfunc'

        if self.id == 1 and self.algo.workers != 1:
            self.parallel = False

    def __dir__(self):
        return ['id','normalised']

# -------------------------------------------------------------------------- DE
class DE:
    """
    Differential Evolution (DE) from scipy.optimize - ID: 1
    . documentation @ Documentation/scipy.optimize.differential_evolution.pdf
    """
    def __init__(self):
        self.strategy = 'best1bin'
        self.maxiter = 1000
        self.popsize = 10
        self.tol = 1e-8
        self.mutation = 0.5
        self.recombination = 0.7
        self.disp = True
        self.polish = False
        self.init = 'latinhypercube'
        self.atol = 1e-15
        self.workers = 1
        self.updating = 'immediate'

        # DO NOT EDIT
        if self.workers != 1:
            self.updating = 'deferred'

# -------------------------------------------------------------------------- DA
class DA:
    """
    Dual Annealing (DA) from Scipy.optimize - ID: 2
    . documentation @ Documentation/scipy.optimize.dual_annealing.pdf
    """
    def __init__(self):
        self.maxiter = 1000
        self.local_search_options = {}
        self.initial_temp = 5230.0
        self.restart_temp_ratio = 2e-5
        self.visit = 2.62
        self.accept = -5.0
        self.maxfun = 1e7
        self.no_local_search = False
        self.x0 = None

# -------------------------------------------------------------------------- LS
class LS:
    """
    Least-Squares (LS) from Scipy.optimize - ID: 3
    . documentation @ Documentation/scipy.optimize.least_squares.pdf
    """
    def __init__(self):
        self.jac = '2-point'
        self.method = 'lm'
        self.bounds = (-np.inf,np.inf)
        self.ftol = 1.0e-8
        self.xtol = 1.0e-8
        self.gtol = 1.0e-8
        self.x_scale = 'jac'
        self.loss = 'linear'
        self.f_scale = 1.0
        self.diff_step = 1e-3
        self.tr_solver = None
        self.tr_options = {}
        self.jac_sparsity = None
        self.max_nfev = 1
        self.verbose = 0

        self.multistart = False
        if self.multistart:
            self.nsamples = 10			# int
            self.multi = 'LHS'			# ['Random','LHS']

        self.usebounds = True

        # DO NOT EDIT
        if self.usebounds:
            if self.method == 'lm':
                self.transform = True
            else:
                self.transform = False

# ------------------------------------------------------------------------ MINI
class MINI:
    """
    Minimize (MINI) from Scipy.optimize - ID: 4
    . documentation @ Documentation/scipy.optimize.minimize.pdf
    """
    def __init__(self):
        self.method = 'Nelder-Mead'
        self.jac = None
        self.hess = None
        self.hessp = None
        self.bounds = None
        self.tol = None
        self.options = {}
        self.options['maxiter'] = 10000
        self.options['disp'] = True

        if self.method == 'Nelder-Mead':
            self.options['maxfev'] = 100000
            self.options['return_all'] = False
            self.options['xatol'] = 1e-8
            self.options['fatol'] = 1e-8
            self.options['adaptive'] = True
        elif self.method == 'Powell':
            self.options['maxfev'] = 100000
            self.options['return_all'] = False
            self.options['xtol'] = 1e-8
            self.options['ftol'] = 1e-8
            self.options['direc'] = None
        elif self.method == 'SLSQP':
            self.options['ftol'] = 1e-8
            self.options['eps'] = 1e-3
            self.options['finite_diff_rel_step'] = None

        self.usebounds = True

        self.multistart = False
        if self.multistart:
            self.nsamples = 10			# int
            self.multi = 'LHS'			# ['Random','LHS','MC]

        # DO NOT EDIT
        meths = ['Nelder-Mead','CG','BFGS','Newton-CG','COBYLA','dogleg',
                 'trust-ncg','trust-krylov','trust-exact']
        if self.usebounds:
            if self.method in meths:
                self.transform = True
            else:
                self.transform = False