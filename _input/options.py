# IMPORT PACKAGES

# IMPORT MODULES
import _input

# --------------------------------------------------------------------- Options
class Options:
    def __init__(self):
        self.info = _input.Info()
        self.method = _input.Method()
        self.optimisation = _input.Optimisation()
        self.material = _input.Material()

    def __dir__(self):
        return ['info','method','optimisation','material']
