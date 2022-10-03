# IMPORT PACKAGES
import numpy as np

# IMPORT MODULES

# ---------------------------------------------------------------------- Method
class Method:
    def __init__(self):
        self.id = 1
        self.residual = 1

        if self.id == 1:
            self.name = 'femu'
            self.description = 'Finite Element Model Updating'
            self.femu = FEMU()
            self.femu.weight = 1.0

        elif self.id == 2:
            self.name = 'vfm'
            self.description = 'Virtual Fields Method'
            self.vfm = VFM()
            self.vfm.weight = 1.0

        elif self.id == 3:
            self.name = 'femu-vfm'
            self.description = 'Finite Element Model Updating + Virtual Fields Method'
            self.femu = FEMU()
            self.vfm = VFM()

            if np.sum([self.femu.weight,self.vfm.weight]) != 1.0:
                print("Error: sum of method's weight is different from 1.0")
                exit()

        if self.id in [1,3]:
            if self.femu.residual == 0:
                self.residual = 2

    def __dir__(self):
        if self.id == 1:
            return ['id','residual','name','femu']
        elif self.id == 2:
            return ['id','residual','name','vfm']
        elif self.id == 3:
            return ['id','residual','name','femu','vfm']

# ------------------------------------------------------------------------ FEMU
class FEMU:
    def __init__(self):
        self.norm = True
        self.weight = 1.0
        self.residual = 2
        self.umat = 'UMAT.obj'
        self.cpus = 1

    def __dir__(self):
        return ['norm','weight','residual','umat','cpus']

# ------------------------------------------------------------------------- VFM
class VFM:
    def __init__(self):
        self.norm = False
        self.weight = 1.0

    def __dir__(self):
        return ['norm','weight']