# IMPORT PACKAGES
import numpy as np
from os import mkdir,listdir,remove

# IMPORT MODULES

# ------------------------------------------------------------------ npy_save()
def npy_save(opti,options,run):
	path = 'Output/npy/'+str(run)+'/'
	try:
		for f in listdir(path):
			remove(path+f)
	except:
		pass
	
	try:
		mkdir(path)
	except:
		pass

	ntests = len(options.info.name)
	nevals = len(opti.variables.all)
	evals = np.arange(1,nevals+1)
	best = opti.best.evolution

	Var_Evol_All = np.column_stack((evals,opti.variables.all))
	Var_Evol_Best = np.column_stack((best,evals,opti.variables.all))

	np.save(path+'Var_Evol_All.npy',Var_Evol_All)
	np.save(path+'Var_Evol_Best.npy',Var_Evol_Best)

	total = opti.objfunc.all.total
	if ntests == 1:
		data = np.column_stack((evals,total))
	else:
		data = np.column_stack((evals,np.sum(total,axis=1)))
	np.save(path+'Obj_Evol_All.npy',data)
	
	total = opti.objfunc.best.total
	if ntests == 1:
		data = np.column_stack((best,evals,total))
	else:
		data = np.column_stack((best,evals,np.sum(total,axis=1)))
	np.save(path+'Obj_Evol_Best.npy',data)

	for n in range(0,ntests):
		name = options.info.name[n]

		total = opti.objfunc.all.total
		femu = opti.objfunc.all.femu
		vfm = opti.objfunc.all.vfm

		if ntests == 1:
			data1 = np.column_stack((evals,total,femu,vfm))
		else:
			data1 = np.column_stack((evals,total[:,n],femu[:,n],vfm[:,n]))

		np.save(path+'Obj_Evol_All_{:s}.npy'.format(name),data)

		total = opti.objfunc.best.total
		femu = opti.objfunc.best.femu
		vfm = opti.objfunc.best.vfm

		if ntests == 1:
			data1 = np.column_stack((best,evals,total,femu,vfm))
		else:
			data1 = np.column_stack((best,evals,total[:,n],femu[:,n],vfm[:,n]))

		np.save(path+'Obj_Evol_Best_{:s}.npy'.format(name),data)
