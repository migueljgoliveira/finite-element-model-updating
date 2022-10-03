# IMPORT MODULES --------------------------------------------------------------
import numpy as np

import _fortran

def vfm(x,opti,options,test,flag):

    # COPY VFM (F2PY) ARGUMENTS
    n = len(x)
    incmax = len(test.data.time)
    nelems = len(test.model.elements)
    ntests = len(options.info.name)
    name = test.info.name
    elprops = options.material.elastic
    thickness = test.model.thickness
    length = test.model.length
    norm = options.method.vfm.norm
    residual = options.method.residual
    expforce = test.data.force
    expstrain = test.data.strain
    maxexpforce = np.max(test.data.force)
    maxexpstrain = np.max(test.data.strain)
    area = test.data.area
    detdfgrd = test.data.detdfgrd
    dfgrd = test.data.dfgrd
    invdfgrd = test.data.invdfgrd
    rot = test.data.rot

    # COMPUTE NUMBER OF RESIDUALS
    m = incmax

    # CALL VFM (F2PY)
    fvec,objfunc = _fortran.vfm(n,m,incmax,nelems,ntests,x,name,elprops,thickness,length,norm,residual,expforce,expstrain,maxexpforce,maxexpstrain,area,detdfgrd,dfgrd,invdfgrd,rot,flag)

    # NORMALIZE BY VFM WEIGHT
    weight = options.method.vfm.weight
    
    fvec = fvec*weight
    objfunc = objfunc*weight

    return fvec,objfunc