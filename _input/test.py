# IMPORT PACKAGES
import numpy as np

# IMPORT MODULES

# ------------------------------------------------------------------------ Test
class Test:
    def __init__(self,name,weight):
        path = 'Input/'
        file1 = '{:s}/{:s}/Model.inp'.format(path,name)
        file2 = '{:s}/{:s}/ExpData.inp'.format(path,name)
        with open(file1,'r') as f:
            fdata = f.readlines()

        self.info  = Info(name,weight)
        self.model = Model(fdata,file1)
        self.data  = Data(fdata,file1,len(self.model.coordinates),file2)

    def __dir__(self):
        return ['info','model','data']

# ------------------------------------------------------------------------ Info
class Info:
    def __init__(self,name,weight):
        self.name  = name
        self.weight  = weight

# ----------------------------------------------------------------------- Model
class Model:
    def __init__(self,fdata,file1):
        self.thickness = None
        self.length = None
        self.orientation = None
        self.coordinates = None
        self.elements = None

        n = 0
        for line in fdata:
            if line.startswith('*Thickness'):
                self.thickness = np.float(fdata[n+1].rstrip('\n'))

            if line.startswith('*Length'):
                self.length = np.float(fdata[n+1].rstrip('\n'))

            if line.startswith('*Orientation'):
                self.orientation = np.float(fdata[n+1].rstrip('\n'))

            elif line.startswith('*Nodes'):
                self.coordinates = np.loadtxt(file1,delimiter=',',skiprows=n+2,
                    max_rows=int(fdata[n+1]))

            elif line.startswith('*Elements'):
                self.elements = np.loadtxt(file1,delimiter=',',skiprows=n+3,
                    max_rows=int(fdata[n+1]))

            n+=1

    def __dir__(self):
        return ['thickness','length','orientation','coordinates','elements']


# ------------------------------------------------------------------------ Data
class Data:
    def __init__(self,fdata,file1,nnodes,file2):
        self.time = None
        self.force = None
        self.displacements = None
        self.strain = None
        self.dfgrad = None
        self.invdfgrd = None
        self.detdfgrd = None
        self.rot = None
        self.area = None

        n = 0
        for line in fdata:
            if line.startswith('*Time'):
                self.time = np.loadtxt(file1,delimiter=',',skiprows=n+2,
                    max_rows=int(fdata[n+1]))
            elif line.startswith('*Reaction'):
                self.force = np.loadtxt(file1,delimiter=',',skiprows=n+2,
                    max_rows=int(fdata[n+1]))

            n+=1

        self.displacements = np.zeros((len(self.time),nnodes,2))

        for i in range(len(self.time)):
            self.displacements[i,:,0] = np.loadtxt(file2,delimiter=',',
                skiprows=1+i*nnodes,max_rows=nnodes,usecols=-2)
            self.displacements[i,:,1] = np.loadtxt(file2,delimiter=',',
                skiprows=1+i*nnodes,max_rows=nnodes,usecols=-1)

    def allocate(self,strain,dfgrd,indfgrd,detdfgrd,rot,area):
        self.strain = strain
        self.dfgrd = dfgrd
        self.invdfgrd = indfgrd
        self.detdfgrd = detdfgrd
        self.rot = rot
        self.area = area

    def __dir__(self):
        return ['time','force','displacements','strain','dfgrd','invdfgrd',
                'detdfgrd','rot','area']

