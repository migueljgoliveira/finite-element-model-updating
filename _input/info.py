# IMPORT PACKAGES
import numpy as np

# IMPORT MODULES

# ------------------------------------------------------------------------ Info
class Info:
    def __init__(self):
        self.calibration = True
        self.name = ['CurvedBar']
        self.weight = [1.0]

        # if np.sum(self.weight) != 1.0:
        # 	print("Error: sum of tests' weight is different from 1.0")
        # 	exit()

    def __dir__(self):
        return ['name','weight']