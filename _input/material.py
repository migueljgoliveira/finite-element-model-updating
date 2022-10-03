# IMPORT PACKAGES
import numpy as np

# IMPORT MODULES

# -------------------------------------------------------------------- Material
class Material:
    def __init__(self):
        self.elastic = [210.0e3, 0.3] # does not work yet
        self.yieldf = 'vonMises'      # does not work yet
        self.hardening = 'swift'      # does not work yet
        self.parameters = {'initial': np.array([ 300.0, 160.00, 0.15]),
                           'lower'  : np.array([ 160.0, 100.00, 0.05]),
                           'upper'  : np.array([ 960.0, 250.00, 0.40])}

    def __dir__(self):
        return ['elastic','yieldf','hardening','parameters']